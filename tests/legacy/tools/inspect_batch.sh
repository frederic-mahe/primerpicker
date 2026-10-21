#!/bin/sh
# Record the Primer Analysis, Hairpins and Dimers screens of many primers
# with inspect_primer.sh, several emulators at a time.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe inspect_batch.sh LIST RUNDIR PAUSE PARALLEL
#
# LIST: one primer per line (a tab and a comment may follow). Each run
# keeps RUNDIR/PRIMER/screens.stripped.txt (screens.txt without the
# splash screens, see strip_splash.py) and screens.txt; the video
# capture, the decoded frames and the copy of PRIMER.EXE are deleted.
# RUNDIR/status.txt gets "PRIMER<TAB>ok" or "PRIMER<TAB>incomplete" per
# run (rerun the incomplete ones with a larger PAUSE). Use PAUSE 1.5 to
# 2.5 and at most 6 parallel runs; primers longer than about 35 bases
# lose keys more often. Run it as a tracked background task.
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
list=$1
rundir=$2
pause=$3
parallel=$4
mkdir -p "$rundir"
export PRIMER_EXE PAUSE="$pause" tools rundir
cut -f1 "$list" | grep -v '^#' | grep . | xargs -P "$parallel" -I{} sh -c '
    sh "$tools/inspect_primer.sh" "$rundir" "$1" >> "$rundir/status.txt" 2>&1
    python3 -I "$tools/strip_splash.py" "$rundir/$1/screens.txt" > "$rundir/$1/screens.stripped.txt"
    rm -rf "$rundir/$1/cap" "$rundir/$1/frames" "$rundir/$1/c"
' sh {}
