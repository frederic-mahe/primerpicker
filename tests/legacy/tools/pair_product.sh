#!/bin/sh
# Record the primer pair and amplified product values that PRIMER.exe
# shows for two primers taken from a molecule.
#
# usage: PRIMER_EXE=/path/to/PRIMER.exe pair_product.sh WORKDIR MOLECULE POS_A POS_B
#
# Primer A (20 bases, the default length) starts at POS_A on the
# molecule; primer B is on the complement strand, at POS_B. Prints
# "MOLECULE<TAB>POS_A<TAB>POS_B<TAB>primer A<TAB>primer B<TAB>Tm A<TAB>Tm B
# <TAB>product size<TAB>product Tm<TAB>annealing temperature"; values
# are empty if keystrokes were lost (retry with a larger PAUSE).
#
# The molecule file must be listed after PRIMER.EXE in the program's
# file list (alphabetical order), as it is selected with two downs.
#
# Keys: space (splash screen), enter (Primer/Create -> from Molecule),
# enter (Load File), down down enter (the molecule), enter (Bp
# position), POS_A, enter, then enter three times (strand N, length
# 20, no extension), s (Save+Exit); enter (from Molecule: the loaded
# molecule is reused), enter (Bp position), POS_B, enter, y (complement
# strand), enter three times, s; right right down enter (Analyze ->
# Primer Pair), 1 enter 2 enter (primers A and B), a (AnalProd).
# The program is slow after s: PAUSE (seconds between keys) must be
# 1.5 or more.
set -eu

tools=$(cd "$(dirname "$0")" && pwd)
workdir=$1
molecule=$2
position_a=$3
position_b=$4
pause=${PAUSE:-1.5}

digits() { printf '%s' "$1" | sed 's/./& /g'; }

name=$(basename "$molecule")
run="$workdir/${name}_${position_a}_${position_b}"
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
autotype -w 4 -p $pause space enter enter down down enter enter $(digits "$position_a") enter enter enter enter s enter enter $(digits "$position_b") enter y enter enter enter s right right down enter 1 enter 2 enter a
dx-capture /V primer.exe
exit
EOC

seconds=$(python3 -c "print(int(60 + 50 * $pause))")
(cd "$run" && SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL "$seconds" dosbox-x -conf run.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" "$run/cap/primer_000.avi" "$run/frames" > /dev/null
python3 -I "$tools/screen_text.py" "$run/frames" > "$run/screens.txt"

screens="$run/screens.txt"
value() { grep -m1 "$1" "$screens" | sed "$2" || true; }
primer_a=$(value "Primer A (bp" "s/.*5'-//; s/-3'.*//")
primer_b=$(value "Primer B (bp" "s/.*5'-//; s/-3'.*//")
tms=$(value "Tm°C  " "s/^│ *Tm°C *//; s/ *│$//; s/  */	/")
size=$(value "Product size:" "s/.*: *//; s/ bps.*//")
product_tm=$(value "Melting temperature:" "s/.*: *//; s/°C.*//")
annealing_tm=$(value "Annealing temperature:" "s/.*: *//; s/°C.*//")
[ -n "$tms" ] || tms="	"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$position_a" "$position_b" \
    "$primer_a" "$primer_b" "$tms" "$size" "$product_tm" "$annealing_tm"
