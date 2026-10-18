#!/bin/sh
# Record the Find Primer Pair list for a molecule.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe find_pair.sh WORKDIR MOLECULE TARGET_A TARGET_B [WEIGHTS]
#   WEIGHTS: five digits for GC mismatch, Tm mismatch, 3' dimers, any
#   dimers and false priming (0 to 5), for instance 05001; default:
#   the built-in weights (all 1)
#
# The molecule file must be listed after PRIMER.EXE. Keys: space, down
# down down enter (Find Pair), enter (Load File), down down enter (the
# molecule), enter (Bp position), TARGET_A, enter x3 (search start,
# extension, length), TARGET_B, enter; with weights: f4, then the five
# digits (one-digit fields move on by themselves); enter x5 (search
# end, extension, length, GC and Tm tolerances) starts the search; a
# key shows the list, and downs scroll it. The list is in screens.txt;
# the counts of possible primers and pairs are on the screen before it.
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
molecule=$2
target_a=$3
target_b=$4
weights=${5:-}
pause=${PAUSE:-2}

digits() { printf '%s' "$1" | sed 's/./& /g'; }
weight_keys=""
[ -n "$weights" ] && weight_keys="f4 $(digits "$weights")"
downs=$(printf 'down %.0s' $(seq 1 12))

name=$(basename "$molecule")
run="$workdir/${name}_${target_a}_${target_b}${weights:+_$weights}"
rm -rf "$run"
mkdir -p "$run/c"
cp "$PRIMER_EXE" "$run/c/PRIMER.EXE"
cp "$molecule" "$run/c/$name"

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
autotype -w 4 -p $pause space down down down enter enter down down enter enter $(digits "$target_a") enter enter enter enter $(digits "$target_b") enter $weight_keys enter enter enter enter enter $downs
dx-capture /V primer.exe
exit
EOC

seconds=$(python3 -c "print(int(150 + 60 * $pause * 2))")
(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$seconds" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
if grep -q 'FIND PRIMER PAIR' "$run/screens.txt"; then
    printf '%s\tok\n' "$run"
else
    printf '%s\tincomplete\n' "$run"
fi
