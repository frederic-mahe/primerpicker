#!/bin/sh
# Record the Mutagenize Primer screen of a primer made from a molecule.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe mutagenize.sh WORKDIR MOLECULE POSITION [KEYS]
#   POSITION: the primer's 5' end (top strand, 20 bases, no extension)
#   KEYS: autotype buttons typed on the Mutagenize screen (',' = one
#   more pause), for instance "a right right t" (overwrite bases) or
#   "f2 down down" (the codon table); the run is named after them
#
# The molecule file must be listed after PRIMER.EXE. Keys: space, enter
# (from Molecule), enter (Load File), down down enter (the molecule),
# two pauses (keys are lost while the molecule map is drawn), enter (Bp
# position), POSITION, enter, enter x3 (strand, length, extension), s
# (Save+Exit), right (Modify), down x5, enter (Mutagenize), then KEYS.
# The screens are in WORKDIR/RUN/screens.txt; prints the run directory
# and "ok" or "incomplete".
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
molecule=$2
position=$3
keys=${4:-}
pause=${PAUSE:-3}

digits() { printf '%s' "$1" | sed 's/./& /g'; }
name=$(basename "$molecule")
tag=$(printf '%s' "$keys" | tr -d ', ')
run="$workdir/${name}_$position${tag:+_$tag}"
rm -rf "$run"
mkdir -p "$run/c"
cp "$PRIMER_EXE" "$run/c/PRIMER.EXE"
cp "$molecule" "$run/c/$name"
line="autotype -w 4 -p $pause space enter enter down down enter , , enter $(digits "$position") enter enter enter enter s , right down down down down down enter , , $keys , ,"
if [ ${#line} -gt 250 ]; then
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

seconds=$(python3 -c "print(int(90 + ${#line} / 4 * $pause * 1.5))")
(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$seconds" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
rm -rf "$run/frames"
if grep -q 'MUTAGENIZE PRIMER' "$run/screens.txt"; then
    printf '%s\tok\n' "$run"
else
    printf '%s\tincomplete\n' "$run"
fi
