#!/bin/sh
# Record the Tm that PRIMER.exe shows for one primer typed with
# "Primer/Create -> Direct Entry".
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe tm_direct_entry.sh WORKDIR PRIMER [EXTENSION]
#
# Prints "EXTENSION<TAB>PRIMER<TAB>result line" (EXTENSION is "-" when
# there is none) and keeps the recording and its decoded screens in
# WORKDIR/EXTENSION_PRIMER/. Runs in parallel are fine, for instance:
#   xargs -P 10 -n1 tm_direct_entry.sh work < primers.txt
#
# Keys: space dismisses the splash screen (the Primer/Create menu is
# then open), down and enter select Direct Entry, the extension is
# typed and confirmed with enter, then the primer. PAUSE is the delay
# between keys, in seconds. If the result is empty, keystrokes were
# lost (the 15-key BIOS buffer overflows on a busy machine): retry
# with PAUSE=1.
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
primer=$2
extension=${3:-}
pause=${PAUSE:-0.3}

keys() { printf '%s' "$1" | tr 'A-Z' 'a-z' | sed 's/./& /g'; }

run="$workdir/${extension:--}_$primer"
rm -rf "$run"
mkdir -p "$run/c"
cp "$PRIMER_EXE" "$run/c/PRIMER.EXE"
length=$((${#extension} + ${#primer}))
seconds=$(python3 -c "print(int(20 + ($length + 6) * $pause * 4))")

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
autotype -w 4 -p $pause space down enter $(keys "$extension") enter $(keys "$primer") enter
dx-capture /V primer.exe
exit
EOC

# DOSBox-X does not exit by itself, and needs SIGKILL
(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$seconds" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
result=$(grep -m1 'Length:.*Tm' "$run/screens.txt" | sed 's/^│ *//; s/ *│$//' || true)
printf '%s\t%s\t%s\n' "${extension:--}" "$primer" "$result"
