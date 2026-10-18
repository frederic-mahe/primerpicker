#!/bin/sh
# Record the Primer Analysis, Hairpins and Dimers screens that PRIMER.exe
# shows for one primer typed with Primer/Create -> Direct Entry.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe inspect_primer.sh WORKDIR PRIMER [EXTENSION]
#
# Keeps the decoded screens in WORKDIR/PRIMER/screens.txt (or
# WORKDIR/EXTENSION_PRIMER/ with a 5' extension; read them with
# parse_inspect.py), and prints "PRIMER<TAB>ok" or
# "PRIMER<TAB>incomplete" (keystrokes lost: retry with a larger PAUSE).
#
# Keys: as tm_direct_entry.sh, then i (Inspect: Primer Analysis
# screen), h (Hairpins screen), esc, d (Dimers screen), esc.
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
primer=$2
extension=${3:-}
pause=${PAUSE:-1}

keys() { printf '%s' "$1" | tr 'A-Z' 'a-z' | sed 's/./& /g'; }

run="$workdir/${extension:+${extension}_}$primer"
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
autotype -w 4 -p $pause space down enter $(keys "$extension") enter $(keys "$primer") enter i h esc d esc
dx-capture /V primer.exe
exit
EOC

seconds=$(python3 -c "print(int(30 + (${#extension} + ${#primer} + 12) * $pause * 3))")
(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$seconds" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
if grep -q 'PRIMER DIMERS' "$run/screens.txt"; then
    printf '%s\tok\n' "$primer"
else
    printf '%s\tincomplete\n' "$primer"
fi
