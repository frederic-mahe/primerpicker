#!/bin/sh
# Record the Multi Primer Analysis screen of 1 to 4 primers.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe multi_primer.sh WORKDIR SELECTION PRIMER...
#   SELECTION: what is typed in the fields Primer A to D, one character
#   per field: a list number (1 to 4) or '_' for an empty field (Enter),
#   for instance 1234, 12__, or 1_2_ (A = 1, C = 2). All four fields are
#   needed: the screen is shown after the last one
#
# The primers are written as .PD2 files (make_pd2.py) and loaded with
# Primer/Create -> Retrieve from disk: typing them with Direct Entry
# needs too many keys (the autotype line must stay under about 250
# characters, or DOSBox-X stops typing at a random point). Keys: space,
# down x6 (the separator counts), enter (Retrieve from disk), down x1,
# enter (P1.PD2); then, for each other primer, enter (the menu stays on
# Retrieve from disk), down xN, enter (PN.PD2); right right down down
# enter (Analyze -> Multi Primer); the SELECTION (one-digit fields move
# on by themselves). The screens are in WORKDIR/RUN/screens.txt (read
# them with parse_multi.py); prints the run directory and "ok" or
# "incomplete".
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
selection=$2
shift 2
if [ ${#selection} -ne 4 ]; then
    echo "SELECTION needs four characters (A to D)" >&2
    exit 1
fi
pause=${PAUSE:-1.8}

run="$workdir/${selection}_$(printf '%s' "$*" | tr ' ' '_')"
rm -rf "$run"
mkdir -p "$run/c"
cp "$PRIMER_EXE" "$run/c/PRIMER.EXE"

keys="space down down down down down down enter"
number=1
for primer in "$@"; do
    python3 -I "$tools/make_pd2.py" "$primer" > "$run/c/P$number.PD2"
    [ "$number" -gt 1 ] && keys="$keys enter"
    keys="$keys $(printf 'down %.0s' $(seq 1 "$number"))enter"
    number=$((number + 1))
done
fields=$(printf '%s' "$selection" | sed 's/./& /g; s/_/enter/g')
line="autotype -w 4 -p $pause $keys right right down down enter $fields"
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

seconds=$(python3 -c "print(int(60 + ${#line} / 4 * $pause * 1.5))")
(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$seconds" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
rm -rf "$run/frames"
if grep -q 'Primer Balance and Cautions' "$run/screens.txt"; then
    printf '%s\tok\n' "$run"
else
    printf '%s\tincomplete\n' "$run"
fi
