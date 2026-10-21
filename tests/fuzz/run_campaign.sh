#!/bin/sh
# Run a fuzzing campaign: every libFuzzer harness, side by side, on this
# machine (proposal §3.10; not run in CI).
#
# usage: sh tests/fuzz/run_campaign.sh MAX_HOURS [PLATEAU_MINUTES [WORKERS]]
#   from the repository root
#   MAX_HOURS: the campaign stops after this time at the latest
#   PLATEAU_MINUTES: a harness stops earlier when its coverage features
#     ("ft:" in the libFuzzer log) have not grown for this long
#     (default 120); the campaign ends when every harness has stopped
#   WORKERS: fuzzing processes in all (default: number of CPUs - 1),
#     shared between harnesses by the weights below
#   CHECK_SECONDS: how often coverage is checked (default 600)
#
# Each harness runs in libFuzzer's fork mode: it goes on after a crash,
# a timeout or an out-of-memory input (each is saved), and merges new
# inputs into its corpus. Corpora are kept between campaigns in
# build/fuzz/corpus/ (seeded from tests/legacy/ the first time; the
# hand-written seeds below are written at every campaign); logs
# and crashing inputs go to build/fuzz/campaigns/DATE/HARNESS/. A
# crashing input is replayed with: build/fuzz/tests/fuzz/HARNESS FILE
set -eu
export DEBUGINFOD_URLS=""  # symbolize from local files only

max_hours=$1
plateau_minutes=${2:-120}
workers=${3:-$(($(nproc) - 1))}
check_seconds=${CHECK_SECONDS:-600}  # how often coverage is checked

# harness:weight (share of the workers)
harnesses="fuzz_find_pair:6 fuzz_cli:5 fuzz_molecule_file:4 fuzz_genbank:3 fuzz_fasta:2 fuzz_site_search:2 fuzz_multi_primer:2 fuzz_pd2:1 fuzz_pd2_cfg:1 fuzz_homology_report:2"
total_weight=0
for entry in $harnesses; do
    total_weight=$((total_weight + ${entry#*:}))
done

root=$(pwd)
build="$root/build/fuzz"
campaign="$build/campaigns/$(date +%Y%m%d%H%M%S)"
max_seconds=$(awk -v hours="$max_hours" 'BEGIN { printf "%d", hours * 3600 }')

cmake --preset fuzz > /dev/null
cmake --build --preset fuzz -j > /dev/null
mkdir -p "$campaign"

# hand-written seeds, (re)written at every campaign so that seeds added
# later reach existing corpora
fixed_seeds() {
    corpus="$build/corpus/$1"
    case $1 in
        fuzz_site_search)
            printf 'GAATTC\nTTGAATTCGG\n' > "$corpus/ecori"
            printf 'GCCNNNNNGGC\nTTGCCATGCAGGCTT\n' > "$corpus/bgli"
            printf 'GGATC\nGGATCGATCC\n' > "$corpus/alwi" ;;
        fuzz_cli)
            printf 'find-pair\n/m\n150\n400' > "$corpus/find_pair"
            printf 'find-pair\n--weights\n05001\n/m\n150\n400' > "$corpus/find_pair_weights"
            printf 'find-enzymes\n--almost\n/m\nAluI' > "$corpus/find_enzymes"
            printf 'analyze\n--extension\nGAATTC\nGTAAAACGACGGCCAGT' > "$corpus/analyze"
            printf 'pair\nCCCTGTCATACCAATCTACCCC\nTGTTATGCGCGTGGGG' > "$corpus/pair"
            printf 'multi\nCTCAAAGGCTCAGACATCAG\nGATCCATGCAATGCTAGCTA\nACGTTGCAGGTCAAGCTTGC\nTTGGCCAAGTCGATCGGATC' \
                > "$corpus/multi"
            printf 'multi\nCTCGATCCGTTCCTAAT\n-\nCGGGGCTCATTGCACTC' > "$corpus/multi_empty_field"
            printf -- '--config\n\\cfg\nfind-pair\n--offsets\n40,40\n/m\n150\n400' > "$corpus/config_find_pair"
            printf 'homology\n//m\n51' > "$corpus/homology"
            printf 'homology\n--false-priming\n--report\nr\n///m\n51' > "$corpus/homology_report"
            printf 'analyze\n--report\nr\nGTAAAACGACGGCCAGT' > "$corpus/analyze_report"
            printf 'pair\n--report\nr\nAAAAAAAAAAAA\nCCCCCCCCCCCC' > "$corpus/pair_report"
            printf 'multi\n--report\nr\nCTCGATCCGTTCCTAAT\n-\nCGGGGCTCATTGCACTC' > "$corpus/multi_report"
            printf 'find-pair\n--report\nr\n/m\n150\n400' > "$corpus/find_pair_report"
            printf 'match\n/m\n51\n400' > "$corpus/match"
            printf 'match\n--report\nr\n/m\n51\n400' > "$corpus/match_report"
            printf 'match\n--top-strand\n/m\n50\n600' > "$corpus/match_top_strand"
            printf 'explain\nTa' > "$corpus/explain"
            # an argument cut by a NUL byte (fuzz_cli.cpp, found 2026-10-20)
            printf 'analyze\n--report\000x\n/\nGGGTCGGG' > "$corpus/nul_report" ;;
        fuzz_multi_primer)
            # fields A to D (0x0F), each 20 bases: the 4-primer recording
            # of tests/legacy/multi/
            python3 -I - "$corpus/multi_1234" <<'EOP'
import sys
primers = ['CTCAAAGGCTCAGACATCAG', 'GATCCATGCAATGCTAGCTA', 'ACGTTGCAGGTCAAGCTTGC', 'TTGGCCAAGTCGATCGGATC']
data = bytes([0x0F])
for primer in primers:
    data += bytes([len(primer) - 2]) + bytes('ACGT'.index(base) for base in primer)
open(sys.argv[1], 'wb').write(data)
EOP
            ;;
        *) ;;
    esac
}

seed() {
    corpus="$build/corpus/$1"
    if [ -d "$corpus" ]; then
        fixed_seeds "$1"
        return
    fi
    mkdir -p "$corpus"
    fixed_seeds "$1"
    case $1 in
        fuzz_pd2) find "$root/tests/legacy" -iname '*.pd2' -exec cp {} "$corpus/" \; ;;
        fuzz_pd2_cfg) find "$root/tests/legacy" -iname '*.cfg' -exec cp {} "$corpus/" \; ;;
        fuzz_site_search | fuzz_cli | fuzz_multi_primer) ;;
        fuzz_find_pair)
            # settings (targets 150 and 400, built-in defaults), then
            # TEST600.GB, one byte per base
            python3 -I - "$root/tests/legacy/pairs/TEST600.GB" "$corpus/test600" <<'EOP'
import re, sys
text = open(sys.argv[1]).read()
sequence = ''.join(re.findall('[acgt]+', text.split('ORIGIN')[1])).upper()
settings = bytes([150, 0, 100, 144, 1, 100, 12, 5, 10, 1, 1, 1, 1, 1])
open(sys.argv[2], 'wb').write(settings + bytes('ACGT'.index(base) for base in sequence))
EOP
            ;;
        *) find "$root/tests/legacy" \( -iname '*.gb' -o -iname '*.fa' -o -iname '*.seq' -o -iname '*.emb' \) \
               -exec cp {} "$corpus/" \; ;;
    esac
}

# a harness-specific option: a dictionary, or the longest useful input
options() {
    case $1 in
        fuzz_genbank | fuzz_fasta | fuzz_molecule_file) printf '%s' "-dict=$root/tests/fuzz/genbank.dict" ;;
        fuzz_cli) printf '%s' "-dict=$root/tests/fuzz/cli.dict" ;;
        fuzz_find_pair) printf '%s' "-max_len=1214" ;;  # 14 bytes of settings, 1200 bases
        fuzz_multi_primer) printf '%s' "-max_len=249" ;;  # 1 + 4 x (1 + 60)
        fuzz_pd2_cfg) printf '%s' "-max_len=128" ;;  # the files have 55 bytes
        *) printf '%s' "-max_len=65536" ;;
    esac
}

# the coverage features of the last status line of a harness
features() {
    sed -n 's/.* ft: \([0-9]*\) .*/\1/p' "$campaign/$1/libfuzzer.log" 2>/dev/null | tail -n 1
}

for entry in $harnesses; do
    harness=${entry%:*}
    jobs=$(( (workers * ${entry#*:} + total_weight - 1) / total_weight ))
    seed "$harness"
    mkdir -p "$campaign/$harness"
    (cd "$campaign/$harness" && exec "$build/tests/fuzz/$harness" "$(options "$harness")" \
        -fork="$jobs" -ignore_crashes=1 -ignore_timeouts=1 -ignore_ooms=1 \
        -max_total_time="$max_seconds" -rss_limit_mb=2048 -timeout=10 \
        -artifact_prefix="$campaign/$harness/" "$build/corpus/$harness" \
        > "$campaign/$harness/libfuzzer.log" 2>&1) &
    echo $! > "$campaign/$harness/pid"
    echo "$(date +%H:%M) $harness: $jobs workers" | tee -a "$campaign/campaign.log"
done

# stop each harness when its coverage has reached a plateau
elapsed=0
while [ "$elapsed" -lt "$max_seconds" ]; do
    sleep "$check_seconds"
    elapsed=$((elapsed + check_seconds))
    running=0
    for entry in $harnesses; do
        harness=${entry%:*}
        pid=$(cat "$campaign/$harness/pid")
        kill -0 "$pid" 2> /dev/null || continue
        running=$((running + 1))
        now=$(features "$harness")
        echo "$(date +%H:%M) $harness ft ${now:-?}" >> "$campaign/$harness/features.log"
        last_change_file="$campaign/$harness/last_change"
        if [ ! -f "$last_change_file" ] || [ "$(cut -d' ' -f1 "$last_change_file")" != "${now:-?}" ]; then
            echo "${now:-?} $elapsed" > "$last_change_file"
        elif [ $((elapsed - $(cut -d' ' -f2 "$last_change_file"))) -ge $((plateau_minutes * 60)) ]; then
            echo "$(date +%H:%M) $harness: no new coverage for $plateau_minutes min, stopped" | \
                tee -a "$campaign/campaign.log"
            # TERM, not INT: a background job of a non-interactive shell
            # ignores SIGINT, and so does libFuzzer started from it
            kill -TERM "$pid" 2> /dev/null || true
        fi
    done
    [ "$running" -eq 0 ] && break
done
wait

echo "campaign: $campaign" | tee -a "$campaign/campaign.log"
for entry in $harnesses; do
    harness=${entry%:*}
    found=$(find "$campaign/$harness" -maxdepth 1 \( -name 'crash-*' -o -name 'leak-*' -o -name 'timeout-*' -o -name 'oom-*' \) | wc -l)
    corpus=$(find "$build/corpus/$harness" -type f | wc -l)
    printf '%-20s %3d inputs found (crash, leak, timeout, oom), ft %s, corpus of %d inputs\n' \
        "$harness" "$found" "$(features "$harness")" "$corpus" | tee -a "$campaign/campaign.log"
done
