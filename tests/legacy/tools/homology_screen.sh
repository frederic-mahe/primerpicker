#!/bin/sh
# Record the Uniqueness or False Priming screen of a primer made from a
# molecule (read the result with parse_homology.py).
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe homology_screen.sh WORKDIR MOLECULE POSITION STRAND KIND
#   STRAND: top or complement; KIND: uniqueness or false_priming
#   CFG: a settings file copied as PD2.CFG (default: none, built-in
#   settings); it is listed before the molecule, one more down
#
# The molecule file must be listed after PRIMER.EXE (alphabetical
# order). Keys: space, enter (from Molecule), enter (Load File), down
# down enter (the molecule), enter (Bp position), POSITION, enter, then
# for the top strand: enter enter enter (strand N, length 20, no
# extension), i (Inspect), u or f; for the complement strand: y enter
# enter enter, s (Save+Exit: this flow returns to the main menu), right
# right enter (Analyze -> Active Primer), u or f.
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
molecule=$2
position=$3
strand=$4
kind=$5
pause=${PAUSE:-2}

digits() { printf '%s' "$1" | sed 's/./& /g'; }
case "$kind" in
    uniqueness) key=u ;;
    false_priming) key=f ;;
    *) echo "KIND must be uniqueness or false_priming" >&2; exit 1 ;;
esac
case "$strand" in
    top) after_position="enter enter enter enter i $key" ;;
    complement) after_position="enter y enter enter enter s right right enter $key" ;;
    *) echo "STRAND must be top or complement" >&2; exit 1 ;;
esac

name=$(basename "$molecule")
cfg=${CFG:-}
run="$workdir/${name}_${position}_${strand}_$kind${cfg:+_$(basename "$cfg" .CFG)}"
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
autotype -w 4 -p $pause space enter enter $downs enter enter $(digits "$position") $after_position
dx-capture /V primer.exe
exit
EOC

seconds=$(python3 -c "print(int(120 + 30 * $pause * 2))")
(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$seconds" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
if grep -q 'BpPos\|is unique' "$run/screens.txt"; then
    printf '%s\tok\n' "$run"
else
    printf '%s\tincomplete\n' "$run"
fi
