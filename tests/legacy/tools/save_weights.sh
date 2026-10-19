#!/bin/sh
# Set the five Find Pair weights in Modify -> System Defaults ->
# Weights, quit and save the defaults: PRIMER.exe writes PD2.CFG.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe save_weights.sh WORKDIR WEIGHTS
#   WEIGHTS: five digits (0 to 5) for GC mismatch, Tm mismatch, 3'
#   dimers, any dimers and false priming, in the order of the form
#
# Keys: space, right (Modify), down x7 (the separator counts), enter
# (System Defaults), down down enter (Weights), the five digits (each
# field moves on by itself; the fifth completes the form and the
# System Defaults menu returns), esc (back to the Modify menu), right
# x3 (Quit), enter, y ("Exit program"), y ("System Defaults Changed But
# Not Saved! SAVE CURRENT DEFAULTS (Y/N) ?"). Prints the path of the
# PD2.CFG written, or "incomplete".
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
weights=$2
pause=${PAUSE:-2}

run="$workdir/weights_$weights"
rm -rf "$run"
mkdir -p "$run/c"
cp "$PRIMER_EXE" "$run/c/PRIMER.EXE"

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
autotype -w 4 -p $pause space right down down down down down down down enter down down enter $(printf '%s' "$weights" | sed 's/./& /g') , esc , right right right , enter , , y , , y
dx-capture /V primer.exe
exit
EOC

(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL 150 dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
rm -rf "$run/frames"
if [ -f "$run/c/PD2.CFG" ]; then
    printf '%s\n' "$run/c/PD2.CFG"
else
    printf '%s\tincomplete\n' "$run"
fi
