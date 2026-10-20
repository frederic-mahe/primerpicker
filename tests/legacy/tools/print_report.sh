#!/bin/sh
# Print a report of PRIMER.exe to disk (P, then Print Control -> Print
# to Disk, file name REPORT) and keep it as WORKDIR/RUN/REPORT.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe print_report.sh WORKDIR KIND ARGUMENTS...
#   inspect PRIMER [EXTENSION]        Primer Analysis of a Direct Entry primer
#                                     (with an optional 5' extension)
#   pair PRIMER_A PRIMER_B            Primer Pair Analysis (Direct Entry)
#   pair_files PRIMER_A PRIMER_B      the same, the primers loaded as .PD2
#                                     files (Retrieve from disk, as
#                                     multi_primer.sh): for primers too long
#                                     to type within one autotype line
#   find_pair MOLECULE TARGET_A TARGET_B
#                                     the Find Primer Pair list
#   homology MOLECULE POSITION KIND   Uniqueness (u) or False Priming (f)
#                                     of a top-strand primer at POSITION
#
# The keys reaching each screen are those of inspect_primer.sh,
# pair_direct.sh, find_pair.sh and homology_screen.sh; then p (the Find
# Pair list first asks "N Pairs found.  Print best: [N]": enter; the
# Uniqueness and False Priming screens have a menu bar: right, enter
# opens Print), down down enter (Print Control -> Print to Disk), the
# file name, enter. CFG: a settings file copied as PD2.CFG (default:
# none, built-in settings; it is listed before the molecule file, which
# takes one more down). A molecule file must be listed after PRIMER.EXE. Prints the run directory and "ok" (REPORT
# written) or "incomplete".
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
kind=$2
shift 2
pause=${PAUSE:-2}

keys() { printf '%s' "$1" | tr 'A-Z' 'a-z' | sed 's/./& /g'; }
to_disk="down down enter , , r e p o r t enter , , , ,"
molecule=""
pd2_files=""
case "$kind" in
    inspect)
        path="space down enter $(keys "${2:-}") enter $(keys "$1") enter i"
        print_keys=", , p , , $to_disk" ;;
    pair)
        path="space down enter enter $(keys "$1") enter s enter enter $(keys "$2") enter s right right down enter 1 enter 2 enter"
        print_keys=", , p , , $to_disk" ;;
    pair_files)
        pd2_files="$1 $2"
        path="space down down down down down down enter down enter enter down down enter right right down enter 1 enter 2 enter"
        print_keys=", , , , p , , $to_disk" ;;
    find_pair)
        molecule=$1
        path="space down down down enter enter down down ${CFG:+down }enter enter $(keys "$2") enter enter enter enter $(keys "$3") enter enter enter enter enter enter , , , , space"
        print_keys=", , p , , enter , , $to_disk" ;;
    homology)
        molecule=$1
        path="space enter enter down down enter , , enter $(keys "$2") enter enter enter enter i , , $3"
        print_keys=", , right , , enter , , $to_disk" ;;
    *)
        echo "KIND must be inspect, pair, pair_files, find_pair or homology" >&2
        exit 1 ;;
esac

cfg=${CFG:-}
run="$workdir/${kind}_$(printf '%s' "$*" | tr ' /' '__' | sed 's/.*tests_legacy_//')${cfg:+_$(basename "$cfg" .CFG)}"
rm -rf "$run"
mkdir -p "$run/c"
cp "$PRIMER_EXE" "$run/c/PRIMER.EXE"
if [ -n "$cfg" ]; then
    cp "$cfg" "$run/c/PD2.CFG"
fi
number=1
for primer in $pd2_files; do
    python3 -I "$tools/make_pd2.py" "$primer" > "$run/c/P$number.PD2"
    number=$((number + 1))
done
if [ -n "$molecule" ]; then
    cp "$molecule" "$run/c/"
fi
# PRINT_KEYS: other keys after the screen is reached (to explore it)
line="autotype -w 4 -p $pause $path ${PRINT_KEYS:-$print_keys}"
if [ ${#line} -gt 240 ]; then
    echo "autotype line too long (${#line} characters)" >&2
    exit 1
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
$line
dx-capture /V primer.exe
exit
EOC

seconds=$(python3 -c "print(int(120 + ${#line} / 4 * $pause * 1.5))")
(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$seconds" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
rm -rf "$run/frames"
if [ -f "$run/c/REPORT" ]; then
    cp "$run/c/REPORT" "$run/REPORT"
    printf '%s\tok\n' "$run"
else
    printf '%s\tincomplete\n' "$run"
fi
