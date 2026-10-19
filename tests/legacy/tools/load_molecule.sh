#!/bin/sh
# Load a molecule file with Primer/Create -> from Molecule -> Load File
# and print what PRIMER.exe says: loaded, too_large ("Sequence too
# large"), insufficient_memory ("Insufficient Memory Error") or unknown.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe load_molecule.sh WORKDIR MOLECULE
#
# The molecule file must be listed after PRIMER.EXE. Keys: space, enter
# (from Molecule), enter (Load File), down down enter (the molecule).
# Prints the molecule, its size in bytes, the size of its sequence text
# (after ORIGIN, or after the FASTA header line) and the result.
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
molecule=$2
pause=${PAUSE:-3}

name=$(basename "$molecule")
run="$workdir/load_$name"
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
autotype -w 4 -p $pause space enter enter down down enter
dx-capture /V primer.exe
exit
EOC

(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL 200 dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
rm -rf "$run/frames"
screens="$run/screens.txt"
if grep -q 'Sequence too large' "$screens"; then
    result=too_large
elif grep -q 'Insufficient Memory' "$screens"; then
    result=insufficient_memory
elif grep -q 'Bp position' "$screens"; then
    result=loaded
else
    result=unknown
fi
section=$(python3 -I - "$molecule" <<'EOP'
import sys
data = open(sys.argv[1], 'rb').read()
if b'ORIGIN' in data:
    start = data.index(b'\n', data.index(b'ORIGIN')) + 1
    end = data.find(b'//', start)
    print(len(data[start:end if end >= 0 else len(data)]))
else:
    print(len(data[data.index(b'\n') + 1:]))
EOP
)
printf '%s\t%s\t%s\t%s\n' "$name" "$(wc -c < "$molecule")" "$section" "$result"
