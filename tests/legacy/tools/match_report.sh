#!/bin/sh
# Print the Find Matching Primer list (Primer/Create -> Match Active
# Primer) to disk, with fewer keys than print_report.sh match (whose
# autotype line is too long for 5-digit positions and targets).
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe match_report.sh WORKDIR MOLECULE POSITION TARGET [CFG]
#   POSITION: 5' end of primer A (top strand, form defaults)
#   TARGET: "Target bp position" of primer B; the other fields of the
#   Define Matching Primer form keep their defaults (those of CFG)
#   CFG: a settings file, copied as PD2.CFG (it is listed before the
#   molecule, which takes one more down)
#   PAUSE: seconds between keys (default 2.5)
#
# Keys: those of print_report.sh match (space, enter enter, the
# molecule, enter, POSITION, enter x4, s, down down enter, TARGET, enter
# x8, p, down down enter, REPORT, enter), with one ',' pause after
# the molecule, after s and after the list is shown. Runs
# tools/probe_keys.sh: the screens are in WORKDIR/RUN/screens.txt, the
# report in WORKDIR/RUN/c/REPORT (copied to WORKDIR/RUN/REPORT). Prints
# the run directory and "ok" or "incomplete" (no REPORT).
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
molecule=$2
position=$3
target=$4
cfg=${5:-}

keys() { printf '%s' "$1" | tr 'A-Z' 'a-z' | sed 's/./& /g'; }
name=$(basename "$molecule" .GB)_${position}_$target${cfg:+_$(basename "$cfg" .CFG)}
files=$molecule
if [ -n "$cfg" ]; then
    mkdir -p "$workdir/cfg_$name"
    cp "$cfg" "$workdir/cfg_$name/PD2.CFG"
    files="$files $workdir/cfg_$name/PD2.CFG"
fi
path="space enter enter down down ${cfg:+down }enter , , enter $(keys "$position")enter enter enter enter s , down down enter $(keys "$target")enter enter enter enter enter enter enter enter , ,"
print="p , down down enter , r e p o r t enter , ,"
# shellcheck disable=SC2086 # files: a list of paths without spaces
"$tools/probe_keys.sh" "$workdir" "$name" "$path $print" $files > /dev/null
run="$workdir/$name"
rm -rf "$workdir/cfg_$name"
if [ -f "$run/c/REPORT" ]; then
    cp "$run/c/REPORT" "$run/REPORT"
    printf '%s\tok\n' "$run"
else
    printf '%s\tincomplete\n' "$run"
fi
