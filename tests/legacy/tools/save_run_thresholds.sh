#!/bin/sh
# Set the two run thresholds in Modify -> System Defaults -> Criteria,
# quit and save the defaults: PRIMER.exe writes PD2.CFG.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe save_run_thresholds.sh WORKDIR RUNS GC_RUNS
#   RUNS: "Reject runs of >= N bases (any)"; GC_RUNS: "Reject runs of
#   >= N G or C at 3'end" (one digit each)
#
# Keys: space, right (Modify), down x7 (the separator counts), enter
# (System Defaults), down enter (Criteria), enter x5 (GC min and max, Tm
# min and max, hairpin cutoff: unchanged), 3 enter 7 enter (the dimer
# thresholds, unchanged: these fields take two digits, so each needs
# Enter), RUNS, GC_RUNS (one-digit fields: each moves on by itself, and
# the last completes the form; the System Defaults menu returns), esc
# (back to the Modify menu), right x3 (Quit), enter, y ("Exit
# program"), y (save the defaults).
# Prints the path of the PD2.CFG written, or "incomplete".
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
runs=$2
gc_runs=$3
pause=${PAUSE:-2}

run="$workdir/runs_${runs}_$gc_runs"
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
autotype -w 4 -p $pause space right down down down down down down down enter down enter enter enter enter enter enter 3 enter 7 enter $runs $gc_runs , esc , right right right , enter , , y , , y
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
