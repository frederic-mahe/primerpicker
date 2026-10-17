"""Decode 80x25 text-mode screenshots (720x400 PNG) into text.

Each 9x16 character cell holds two colours; its pixel mask is looked up
in a glyph table made from `vga_glyphs.png`, a screenshot of the 256
code page 437 characters written from the top-left corner (see
`vga_glyphs.conf`). Cells that match no glyph, such as the cell under
the blinking cursor, decode as '?'.

Usage: python3 -I screen_text.py FRAMES_DIR
prints the distinct screens of a recording, in order (screens that
differ only by '?' cells count as one).
"""
import sys
from pathlib import Path

import numpy as np
from PIL import Image

COLUMNS, ROWS = 80, 25
CELL_WIDTH, CELL_HEIGHT = 9, 16
GLYPHS = Path(__file__).with_name('vga_glyphs.png')
UNKNOWN = '?'


def cells(png_path):
    """Yield (row, column, cell) with cell an array of 24-bit colours."""
    rgb = np.asarray(Image.open(png_path).convert('RGB')).astype(np.int32)
    colours = (rgb[..., 0] << 16) | (rgb[..., 1] << 8) | rgb[..., 2]
    for row in range(ROWS):
        for column in range(COLUMNS):
            yield row, column, colours[row * CELL_HEIGHT:(row + 1) * CELL_HEIGHT,
                                       column * CELL_WIDTH:(column + 1) * CELL_WIDTH]


def load_glyphs(png_path=GLYPHS):
    """Map each glyph's foreground mask to its character."""
    glyphs = {}
    for row, column, cell in cells(png_path):
        code = row * COLUMNS + column
        if code > 255:
            break
        mask = (cell != 0).tobytes()
        char = bytes([code]).decode('cp437')
        # several codes share a glyph (0, 32 and 255 are all blank):
        # keep a printable character
        if mask not in glyphs or (char.isprintable() and not glyphs[mask].isprintable()):
            glyphs[mask] = char
    return glyphs


def decode(png_path, glyphs):
    lines = [[' '] * COLUMNS for _ in range(ROWS)]
    for row, column, cell in cells(png_path):
        colours = np.unique(cell)
        if len(colours) == 1:
            continue
        lines[row][column] = UNKNOWN
        if len(colours) > 2:
            continue
        for colour in colours:
            char = glyphs.get((cell == colour).tobytes(), ' ')
            if char != ' ':
                lines[row][column] = char
                break
    return '\n'.join(''.join(line).rstrip() for line in lines)


def distinct_screens(frames_dir):
    glyphs = load_glyphs()
    seen = set()
    for png in sorted(Path(frames_dir).glob('*.png')):
        text = decode(png, glyphs)
        key = text.replace(UNKNOWN, ' ')
        if key in seen:
            continue
        seen.add(key)
        yield png.name, text


if __name__ == '__main__':
    for name, text in distinct_screens(sys.argv[1]):
        print(f'==== {name}')
        print(text)
