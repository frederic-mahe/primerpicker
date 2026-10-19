#!/bin/sh
# Run Find Primer Pair with chosen search windows, to probe its limits
# (100 possible primers, 1000 pairs).
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe find_pair_limits.sh WORKDIR MOLECULE TARGET_A OFFSET_A TARGET_B OFFSET_B SECONDS
#   OFFSET_A: "Start search at: Target -"; OFFSET_B: "Extend search
#   to Target +"
#   SECONDS: time limit of the emulator (searches capped at 1000 pairs
#   never end: the program stops reading keys)
#   CFG: a settings file copied as PD2.CFG (default: none, built-in
#   criteria)
#   EXTRA: autotype buttons typed after the form, ',' being one more
#   pause (for instance a key 30 s after the search starts)
#   CHECK: seconds after which the run stops unless the search started
#   or a message appeared (default 300)
#   PAUSE: seconds between keys (default 4: keys get lost while the
#   molecule map is drawn)
#
# The molecule file must be listed after PRIMER.EXE (and PD2.CFG).
# Keys: space, down x3, enter (Find Pair), enter (Load File), downs and
# enter (the molecule), enter (Bp position), TARGET_A, enter, OFFSET_A,
# enter, enter (length), TARGET_B, enter, OFFSET_B, enter x4 (extension,
# length, GC and Tm tolerances: the last Enter starts the search). The
# 5' extension field after OFFSET_A is skipped by Enter; an offset of 3
# digits fills its field, which then moves on by itself (no Enter). The
# screens are in screens.txt, the counts and messages are printed.
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
molecule=$2
pause=${PAUSE:-4}

digits() { printf '%s' "$1" | sed 's/./& /g'; }
offset() { if [ ${#1} -ge 3 ]; then digits "$1"; else printf '%s enter' "$(digits "$1")"; fi; }

name=$(basename "$molecule")
cfg=${CFG:-}
run="$workdir/${name}_$3_$4_$5_$6${cfg:+_$(basename "$cfg" .CFG)}"
rm -rf "$run"
mkdir -p "$run/c"
cp "$PRIMER_EXE" "$run/c/PRIMER.EXE"
cp "$molecule" "$run/c/$name"
downs="down down"
if [ -n "$cfg" ]; then
    cp "$cfg" "$run/c/PD2.CFG"
    downs="down down down"
fi

cat > "$run/run.conf" <<EOC
[sdl]
output=surface
[dosbox]
machine=svga_s3
captures=cap
[cpu]
cycles=max
[autoexec]
mount c c
c:
autotype -w 4 -p $pause space down down down enter enter $downs enter enter $(digits "$3") enter $(offset "$4") enter enter $(digits "$5") enter $(offset "$6") enter enter enter enter ${EXTRA:-}
dx-capture /V primer.exe
exit
EOC

(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$7" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) &
emulator=$!
# stop early when keys were lost (the form never completed)
sleep "${CHECK:-300}"
mkdir -p "$run/check"
cp "$run/cap/primer_000.avi" "$run/check/partial.avi"
python3 -I "$tools/zmbv.py" "$run/check/partial.avi" "$run/check/frames" > /dev/null 2>&1 || true
python3 -I "$tools/screen_text.py" "$run/check/frames" > "$run/check/screens.txt" 2>/dev/null || true
if ! grep -q 'SEARCH MATCHING PAIRS\|More than\|No Matching' "$run/check/screens.txt"; then
    pkill -KILL -f "dosbox-x -conf run.conf" || true
fi
rm -rf "$run/check"
wait "$emulator" || true

python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
if grep -q 'SEARCH MATCHING PAIRS\|More than\|No Matching' "$run/screens.txt"; then
    printf '%s\tok\n' "$run"
else
    printf '%s\tincomplete\n' "$run"
fi
grep -o 'Primer [AB]: *[0-9]* possible primers\|Pairs: *[0-9]* possible\|More than[^║▒]*\|No Matching Primers Found' \
    "$run/screens.txt" | sed 's/  */ /g' | sort | uniq -c
