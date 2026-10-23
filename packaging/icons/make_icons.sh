#!/bin/sh
# Makes the icon files of PrimerPicker in this folder, from the two
# drawings of make_icon.py: primerpicker.svg (64 px and larger) and
# primerpicker-small.svg (16 to 48 px). The generated files are kept in
# the repository; run this again after changing make_icon.py.
#
# Needs Python 3 with fontTools, the Liberation Sans font, Inkscape and
# ImageMagick 7.
#
# usage: sh packaging/icons/make_icons.sh
set -eu

here=$(cd "$(dirname "$0")" && pwd)
cd "$here"

python3 make_icon.py .
mkdir -p png
# each size is drawn from the vectors (not scaled down from a large image)
for size in 16 24 32 48 64; do
    inkscape primerpicker-small.svg --export-type=png --export-filename="png/primerpicker-small-$size.png" -w "$size" -h "$size"
done
for size in 64 128 256 512 1024; do
    inkscape primerpicker.svg --export-type=png --export-filename="png/primerpicker-$size.png" -w "$size" -h "$size"
done
# PNG metadata (dates) would change the files at each run
for file in png/*.png; do
    magick "$file" -strip "$file"
done

# Windows: the small drawing up to 48 px, the full one from 64 px
magick png/primerpicker-small-16.png png/primerpicker-small-24.png png/primerpicker-small-32.png \
    png/primerpicker-small-48.png png/primerpicker-64.png png/primerpicker-128.png png/primerpicker-256.png \
    primerpicker.ico
# macOS
python3 make_icns.py primerpicker.icns png
rm png/primerpicker-512.png png/primerpicker-1024.png png/primerpicker-small-64.png
