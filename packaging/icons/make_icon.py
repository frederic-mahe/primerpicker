"""Draws the PrimerPicker icons (chosen 2026-10-23): primerpicker.svg,
"PRIMER", a left half-arrow over a right half-arrow (as the sign ⇋),
"PICKER", for 64 px and larger; primerpicker-small.svg, the arrows
alone, thicker, for smaller sizes (the words are unreadable there).
The words are turned into outlines (Liberation Sans Bold, SIL OFL
1.1), so the SVG files need no font. Run by make_icons.sh.

usage: python3 make_icon.py OUTPUT_DIR
"""
import sys
from pathlib import Path

from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

FONT = "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"
CANVAS = 1024
# macOS grid: an 824 px rounded square centred on the 1024 px canvas
TILE, TILE_RADIUS = 824, 185
LEFT, RIGHT = 182, 842            # content width 660
CAP_HEIGHT = 150                  # at most; less if the words are too wide
SHAFT, GAP = 50, 22               # arrow thickness, gap between shafts
BARB = 120                        # barb height (45 degrees)
SHIFT = 40                        # horizontal offset of the two arrows

WORDS_COLOUR = "#1f4e79"          # navy
ARROWS_COLOUR = "#e8891c"         # orange


def word_path(font, word, cap_height, centre_x, baseline):
    glyphs = font.getGlyphSet()
    cmap = font.getBestCmap()
    scale = cap_height / font["OS/2"].sCapHeight
    names = [cmap[ord(char)] for char in word]
    width = sum(glyphs[name].width for name in names) * scale
    pen = SVGPathPen(glyphs)
    x = centre_x - width / 2
    for name in names:
        # font units, y up -> SVG pixels, y down
        glyphs[name].draw(TransformPen(pen, (scale, 0, 0, -scale, x, baseline)))
        x += glyphs[name].width * scale
    return pen.getCommands(), width


def harpoon(left, right, shaft_bottom, shaft=SHAFT, barb=BARB):
    """Half-arrow pointing left, barb up; its bottom edge is shaft_bottom."""
    width = shaft * 2 ** 0.5
    return [(left, shaft_bottom), (left + barb, shaft_bottom - barb),
            (left + barb + width, shaft_bottom - barb),
            (left + width + shaft, shaft_bottom - shaft),
            (right, shaft_bottom - shaft), (right, shaft_bottom)]


def points(polygon):
    return " ".join(f"{x:.1f},{y:.1f}" for x, y in polygon)


def cap_height_to_fit(font):
    # the wider word spans the content width
    widest = max(word_path(font, word, CAP_HEIGHT, 0, 0)[1] for word in ("PRIMER", "PICKER"))
    return min(CAP_HEIGHT, CAP_HEIGHT * (RIGHT - LEFT) / widest)


def icon(words_colour, arrows_colour, font):
    cap_height = cap_height_to_fit(font)
    arrows_height = 2 * BARB + GAP
    spacing = 56
    total = 2 * cap_height + 2 * spacing + arrows_height
    top = (CANVAS - total) / 2
    primer, _ = word_path(font, "PRIMER", cap_height, CANVAS / 2, top + cap_height)
    picker, _ = word_path(font, "PICKER", cap_height, CANVAS / 2, CANVAS - top)
    upper_bottom = top + cap_height + spacing + BARB
    upper = harpoon(LEFT, RIGHT - SHIFT, upper_bottom)
    # the right half-arrow: the left one turned by 180 degrees around
    # the middle of the gap between the shafts
    centre_y = upper_bottom + GAP / 2
    lower = [(CANVAS - x, 2 * centre_y - y) for x, y in upper]
    corner = (CANVAS - TILE) / 2
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{CANVAS}" height="{CANVAS}" viewBox="0 0 {CANVAS} {CANVAS}">
  <rect x="{corner}" y="{corner}" width="{TILE}" height="{TILE}" rx="{TILE_RADIUS}" fill="#ffffff" stroke="#d5dbe1" stroke-width="4"/>
  <g fill="{words_colour}">
    <path d="{primer}"/>
    <path d="{picker}"/>
  </g>
  <g fill="{arrows_colour}">
    <polygon points="{points(upper)}"/>
    <polygon points="{points(lower)}"/>
  </g>
</svg>
"""


# Small sizes (16 to 48 px): the arrows alone, thicker, on a tile that
# fills the canvas
SMALL_TILE, SMALL_RADIUS = 960, 200


def small_tile():
    corner = (CANVAS - SMALL_TILE) / 2
    return (f'<rect x="{corner}" y="{corner}" width="{SMALL_TILE}" height="{SMALL_TILE}" '
            f'rx="{SMALL_RADIUS}" fill="#ffffff" stroke="#d5dbe1" stroke-width="12"/>')


def arrows_icon(arrows_colour):
    shaft, barb, gap, shift, left = 120, 250, 56, 70, 120
    right = CANVAS - left
    upper_bottom = (CANVAS - gap) / 2
    upper = harpoon(left, right - shift, upper_bottom, shaft, barb)
    lower = [(CANVAS - x, CANVAS - y) for x, y in upper]
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{CANVAS}" height="{CANVAS}" viewBox="0 0 {CANVAS} {CANVAS}">
  {small_tile()}
  <g fill="{arrows_colour}">
    <polygon points="{points(upper)}"/>
    <polygon points="{points(lower)}"/>
  </g>
</svg>
"""


def main():
    output = Path(sys.argv[1])
    output.mkdir(parents=True, exist_ok=True)
    font = TTFont(FONT)
    (output / "primerpicker.svg").write_text(icon(WORDS_COLOUR, ARROWS_COLOUR, font))
    (output / "primerpicker-small.svg").write_text(arrows_icon(ARROWS_COLOUR))


if __name__ == "__main__":
    main()
