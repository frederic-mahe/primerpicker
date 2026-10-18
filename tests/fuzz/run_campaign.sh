#!/bin/sh
# Run a fuzzing campaign: every libFuzzer harness, side by side, for a
# given time, on this machine (proposal §3.10; not run in CI).
#
# usage: sh tests/fuzz/run_campaign.sh HOURS [JOBS_PER_HARNESS]
#   from the repository root; JOBS_PER_HARNESS defaults to 4
#
# Builds the "fuzz" preset, seeds each corpus from tests/legacy/ (and
# keeps corpora between campaigns in build/fuzz/corpus/), and writes
# logs and crashing inputs to build/fuzz/campaigns/DATE/. A crashing
# input is reproduced with: build/fuzz/tests/fuzz/HARNESS FILE
set -eu
export DEBUGINFOD_URLS=""  # symbolize from local files only

hours=$1
jobs=${2:-4}
seconds=$(awk -v hours="$hours" 'BEGIN { printf "%d", hours * 3600 }')
root=$(pwd)
build="$root/build/fuzz"
campaign="$build/campaigns/$(date +%Y%m%d%H%M%S)"
harnesses="fuzz_genbank fuzz_fasta fuzz_molecule_file fuzz_pd2 fuzz_site_search"

cmake --preset fuzz > /dev/null
cmake --build --preset fuzz -j > /dev/null
mkdir -p "$campaign"

seed() {
    corpus="$build/corpus/$1"
    mkdir -p "$corpus"
    case $1 in
        fuzz_pd2) find "$root/tests/legacy" -iname '*.pd2' -exec cp {} "$corpus/" \; ;;
        fuzz_site_search)
            printf 'GAATTC\nTTGAATTCGG\n' > "$corpus/ecori"
            printf 'GCCNNNNNGGC\nTTGCCATGCAGGCTT\n' > "$corpus/bgli"
            printf 'GGATC\nGGATCGATCC\n' > "$corpus/alwi" ;;
        *) find "$root/tests/legacy" \( -iname '*.gb' -o -iname '*.fa' -o -iname '*.seq' -o -iname '*.emb' \) \
               -exec cp {} "$corpus/" \; ;;
    esac
}

for harness in $harnesses; do
    seed "$harness"
    mkdir -p "$campaign/$harness"
    dictionary=""
    case $harness in
        fuzz_genbank | fuzz_fasta | fuzz_molecule_file) dictionary="-dict=$root/tests/fuzz/genbank.dict" ;;
    esac
    (cd "$campaign/$harness" && "$build/tests/fuzz/$harness" ${dictionary:+"$dictionary"} \
        -max_total_time="$seconds" -jobs="$jobs" -workers="$jobs" -max_len=65536 -rss_limit_mb=2048 \
        -timeout=10 -artifact_prefix="$campaign/$harness/" "$build/corpus/$harness" \
        > "$campaign/$harness/libfuzzer.log" 2>&1) &
done
wait

echo "campaign: $campaign"
for harness in $harnesses; do
    crashes=$(find "$campaign/$harness" -maxdepth 1 \( -name 'crash-*' -o -name 'leak-*' -o -name 'timeout-*' -o -name 'oom-*' \) | wc -l)
    corpus=$(find "$build/corpus/$harness" -type f | wc -l)
    printf '%-20s %3d crashing inputs, corpus of %d inputs\n' "$harness" "$crashes" "$corpus"
done
