#!/bin/sh
# Type a key sequence into PRIMER.exe, Ctrl keys included, and record
# the screens with their colours (the picking screen of From Molecule
# and Modify -> Position highlights the primer with a background
# colour).
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe probe_picking.sh WORKDIR NAME KEYS [FILE...]
#   the arguments of probe_keys.sh; in KEYS, f5 to f10 stand for
#   Ctrl+Home, Ctrl+End, Ctrl+PgUp, Ctrl+PgDn, Ctrl+Left and Ctrl+Right
#   (CTRLKEYS.COM, built from ctrlkeys.S with GNU as and ld, is loaded
#   from drive D: before PRIMER.EXE, so that the Load File list of C: is
#   that of probe_keys.sh); PAUSE, SECONDS_MAX as for probe_keys.sh
#
# Writes WORKDIR/NAME/screens.txt (screen_text.py) and
# WORKDIR/NAME/colours.txt (screen_colours.py: every change of screen,
# with the colours of each row). Prints the run directory and the new
# files of c/, as probe_keys.sh.
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
name=$2
keys=$3
shift 3
pause=${PAUSE:-2.5}

run="$workdir/$name"
rm -rf "$run"
mkdir -p "$run/c" "$run/d" "$run/build"
as --32 -o "$run/build/ctrlkeys.o" "$tools/ctrlkeys.S"
ld -m elf_i386 -Ttext 0x100 --oformat binary -o "$run/d/CTRLKEYS.COM" "$run/build/ctrlkeys.o"
rm -rf "$run/build"
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
mount d d
d:\\ctrlkeys
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
python3 -I "$tools/screen_colours.py" "$run/frames" > "$run/colours.txt"
rm -rf "$run/frames"
new=$(ls "$run/c" | grep -vxF "$before" || true)
if [ -n "$new" ]; then
    printf '%s\tnew: %s\n' "$run" "$(printf '%s' "$new" | tr '\n' ' ')"
else
    printf '%s\tno new file\n' "$run"
fi
