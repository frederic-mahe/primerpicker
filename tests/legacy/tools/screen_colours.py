"""Decode 80x25 text-mode screenshots with their colours.

The same decoding as screen_text.py, plus each cell's foreground and
background colour, named after the 16 colours of the VGA text-mode
palette (DOSBox-X's default palette). Cells of a single colour (blank)
have only a background.

Usage: python3 -I screen_colours.py FRAMES_DIR [FIRST_ROW LAST_ROW]
prints the screens of a recording (text and colours), in order,
leaving out a screen identical to the one before it (unlike
screen_text.py, a screen seen earlier is printed again: moving back
to an earlier state shows). For each screen: the text, then one line per row of colour
runs, "row R: C1-C2 FG/BG ...", columns and rows counted from 0
(rows FIRST_ROW to LAST_ROW only, when given).
"""
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # python3 -I leaves it out
import screen_text  # noqa: E402

# the 16 text-mode colours of the VGA palette (6-bit DAC values x 4,
# as DOSBox-X writes them)
VGA_COLOURS = {
    0x000000: 'black', 0x0000aa: 'blue', 0x00aa00: 'green', 0x00aaaa: 'cyan',
    0xaa0000: 'red', 0xaa00aa: 'magenta', 0xaa5500: 'brown', 0xaaaaaa: 'lightgrey',
    0x555555: 'darkgrey', 0x5555ff: 'lightblue', 0x55ff55: 'lightgreen',
    0x55ffff: 'lightcyan', 0xff5555: 'lightred', 0xff55ff: 'lightmagenta',
    0xffff55: 'yellow', 0xffffff: 'white',
}


def name(colour):
    return VGA_COLOURS.get(int(colour), f'#{int(colour):06x}')


def cell_colours(cell, glyphs):
    """(foreground, background) names of a cell; foreground '' if blank."""
    colours = np.unique(cell)
    if len(colours) == 1:
        return '', name(colours[0])
    if len(colours) > 2:
        return '?', '?'
    for index, colour in enumerate(colours):
        if glyphs.get((cell == colour).tobytes(), ' ') != ' ':
            return name(colour), name(colours[1 - index])
    return '?', '?'


def colour_rows(png_path, glyphs):
    rows = [[None] * screen_text.COLUMNS for _ in range(screen_text.ROWS)]
    for row, column, cell in screen_text.cells(png_path):
        rows[row][column] = cell_colours(cell, glyphs)
    return rows


def runs(row):
    """(first, last, foreground, background) of each run; a blank cell
    joins a run of the same background."""
    found = []
    for column, (foreground, background) in enumerate(row):
        if found and found[-1][3] == background and (foreground in ('', found[-1][2]) or found[-1][2] == ''):
            first, _, previous, _ = found[-1]
            found[-1] = (first, column, previous or foreground, background)
            continue
        found.append((column, column, foreground, background))
    return found


def main(frames_dir, first_row=0, last_row=screen_text.ROWS - 1):
    glyphs = screen_text.load_glyphs()
    previous = None
    for png in sorted(Path(frames_dir).glob('*.png')):
        text = screen_text.decode(png, glyphs)
        rows = colour_rows(png, glyphs)
        key = (text.replace(screen_text.UNKNOWN, ' '), str(rows[first_row:last_row + 1]))
        if key == previous:
            continue
        previous = key
        print(f'==== {png.name}')
        print(text)
        for number in range(first_row, last_row + 1):
            described = ' '.join(f'{first}-{last} {fore or "-"}/{back}' for first, last, fore, back in runs(rows[number]))
            print(f'row {number}: {described}')


if __name__ == '__main__':
    arguments = sys.argv[1:]
    if len(arguments) == 3:
        main(arguments[0], int(arguments[1]), int(arguments[2]))
    else:
        main(arguments[0])
