#!/bin/sh
# Record the True Sites or Almost Sites screen of one enzyme.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe find_enzyme.sh WORKDIR MOLECULE true|almost ENZYME
#   ENZYME: an enzyme of the default menu (enzymes/enzyme_table.tsv,
#   default_menu = Y); it is reached with downs from the top of the
#   menu, which lists these enzymes in alphabetical order
#   PAUSE: seconds between keys (default 3: keys are lost while the
#   molecule map is drawn); DOWNS: downs that scroll the result list
#   (default 6)
#
# The molecule file must be listed after PRIMER.EXE. Keys: space, right
# x3, enter (Molecule-Ops -> Find Enzymes), enter (Load File), down down
# enter (the molecule), enter (Almost Sites) or down enter (True Sites),
# the downs to the enzyme, enter. The result screens are in screens.txt.
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
molecule=$2
kind=$3
enzyme=$4
pause=${PAUSE:-3}

index=$(awk -F'\t' -v enzyme="$enzyme" '
    !/^#/ && $5 == "Y" { if ($1 == enzyme) { print n; found = 1; exit } n++ }
    END { if (!found) exit 1 }' "$tools/../enzymes/enzyme_table.tsv") || {
    echo "$enzyme is not in the default menu" >&2
    exit 1
}
case $kind in
    true) kind_keys="down enter" ;;
    almost) kind_keys="enter" ;;
    *) echo "kind: true or almost" >&2; exit 1 ;;
esac
repeat() { [ "$2" -gt 0 ] && printf "$1 %.0s" $(seq 1 "$2"); true; }

name=$(basename "$molecule")
run="$workdir/${name}_${kind}_${enzyme}"
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
autotype -w 4 -p $pause space right right right enter enter down down enter $kind_keys $(repeat down "$index") enter $(repeat down "${DOWNS:-6}")
dx-capture /V primer.exe
exit
EOC

seconds=$(python3 -c "print(int(120 + $pause * (20 + $index + ${DOWNS:-6}) * 1.5))")
(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$seconds" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
if grep -q "Search for: $enzyme " "$run/screens.txt"; then
    printf '%s\tok\n' "$run"
else
    printf '%s\tincomplete\n' "$run"
fi
