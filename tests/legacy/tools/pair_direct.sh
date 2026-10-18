#!/bin/sh
# Record the Primer Pair Analysis screen of two primers typed with
# Primer/Create -> Direct Entry.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe pair_direct.sh WORKDIR PRIMER_A PRIMER_B
#
# Keeps the decoded screens in WORKDIR/A_B/screens.txt (read them with
# parse_pair.py), and prints "A<TAB>B<TAB>ok" or "...incomplete" (keys
# lost: retry with a larger PAUSE; 4 runs in parallel with PAUSE=1.8
# is safe).
#
# Keys: space (splash), down enter (Direct Entry), enter (no
# extension), primer A, enter, s (Save+Exit); enter (the menu stays on
# Direct Entry), enter, primer B, enter, s; right right down enter
# (Analyze -> Primer Pair), 1 enter 2 enter.
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
primer_a=$2
primer_b=$3
pause=${PAUSE:-1.8}

keys() { printf '%s' "$1" | tr 'A-Z' 'a-z' | sed 's/./& /g'; }

run="$workdir/${primer_a}_$primer_b"
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
autotype -w 4 -p $pause space down enter enter $(keys "$primer_a") enter s enter enter $(keys "$primer_b") enter s right right down enter 1 enter 2 enter
dx-capture /V primer.exe
exit
EOC

length=$((${#primer_a} + ${#primer_b}))
seconds=$(python3 -c "print(int(40 + ($length + 25) * $pause * 2.5))")
(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$seconds" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
if grep -q 'Any adj homol bases' "$run/screens.txt"; then
    printf '%s\t%s\tok\n' "$primer_a" "$primer_b"
else
    printf '%s\t%s\tincomplete\n' "$primer_a" "$primer_b"
fi
