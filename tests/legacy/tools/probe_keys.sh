#!/bin/sh
# Type any key sequence into PRIMER.exe and record the screens.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe probe_keys.sh WORKDIR NAME KEYS [FILE...]
#   NAME: the run directory (WORKDIR/NAME)
#   KEYS: autotype buttons, from the splash screen on (its first key is
#   taken by the splash screen), ',' being one more pause; for instance
#   "space down enter enter c a g t enter s"
#   FILE: files copied next to PRIMER.EXE (molecules, .PD2 files; a
#   file named PD2.CFG is the settings file). Molecules must be named to
#   sort after PRIMER.EXE in the Load File list.
#   PAUSE: seconds between keys (default 2.5); SECONDS_MAX: emulator
#   time limit (default from the number of keys)
#
# The screens are in WORKDIR/NAME/screens.txt; every file PRIMER.exe
# writes stays in WORKDIR/NAME/c (for instance a report printed to
# disk, or PD2.CFG). Prints the run directory and the files of c/ that
# were not copied in ("new: ..."), or "no new file".
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
name=$2
keys=$3
shift 3
pause=${PAUSE:-2.5}

run="$workdir/$name"
rm -rf "$run"
mkdir -p "$run/c"
cp "$PRIMER_EXE" "$run/c/PRIMER.EXE"
for file in "$@"; do
    cp "$file" "$run/c/"
done
before=$(ls "$run/c")
line="autotype -w 4 -p $pause $keys"
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

seconds=${SECONDS_MAX:-$(python3 -c "print(int(90 + ${#line} / 4 * $pause * 1.5))")}
(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$seconds" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"
rm -rf "$run/frames"
new=$(ls "$run/c" | grep -vxF "$before" || true)
if [ -n "$new" ]; then
    printf '%s\tnew: %s\n' "$run" "$(printf '%s' "$new" | tr '\n' ' ')"
else
    printf '%s\tno new file\n' "$run"
fi
