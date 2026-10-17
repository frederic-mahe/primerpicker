#!/bin/sh
# Remake vga_glyphs.png (the glyph table of screen_text.py).
# usage: make_vga_glyphs.sh WORKDIR
set -eu
tools=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$1/c"
cd "$1"
rm -rf cap frames
cp "$tools/vga_glyphs.conf" .
# GLYPHS.COM is a 24-byte DOS program that writes characters 0-255
# (grey on black) from the top-left corner of the text screen, then
# waits for a key:
#   mov ax, 0xb800 / mov es, ax / xor di, di    ; es:di = text screen
#   mov ah, 0x07 / xor al, al / mov cx, 256     ; grey on black, from 0
#   stosw / inc al / loop -5                    ; write 256 characters
#   xor ax, ax / int 0x16 / ret                 ; wait for a key, exit
printf '\270\000\270\216\300\061\377\264\007\060\300\271\000\001\253\376\300\342\373\061\300\315\026\303' > c/GLYPHS.COM
(SDL_VIDEODRIVER=dummy SDL_AUDIODRIVER=dummy \
    timeout -s KILL 12 dosbox-x -conf vga_glyphs.conf \
    -nopromptfolder -fastlaunch > /dev/null 2>&1) || true
python3 -I "$tools/zmbv.py" cap/glyphs_000.avi frames > /dev/null
# frame 0 is the DOS prompt; the next one shows the characters
cp "$(ls frames/*.png | sed -n 2p)" "$tools/vga_glyphs.png"
