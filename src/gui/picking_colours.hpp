#pragma once

#include "core/colour/colour.hpp"

namespace primerpicker::gui {

// The colours of the picking view (From Molecule, Modify -> Position),
// in one place, a light and a dark variant. PRIMER.exe shows the
// primer's window with a background colour (light grey when the primer
// fails the criteria, green when it meets them); here the two
// backgrounds are a neutral grey and a blue, which stay apart with
// colour-vision deficiencies, and the window has a second cue, its
// border: dashed when it fails, solid when it meets the criteria (and
// the verdict is written in words). Checked by tests/gui/
// test_picking_colours.cpp: WCAG AA contrast of the bases, 3:1 for the
// borders, CIE76 differences of the backgrounds with simulated colour
// vision.
struct PickingPalette {
    colour::Srgb view_background;
    colour::Srgb text;  // the bases, on the page and in the window
    colour::Srgb fails;
    colour::Srgb meets;
    colour::Srgb fails_border;
    colour::Srgb meets_border;
};

constexpr auto light_picking_palette = PickingPalette{
    .view_background = {.red = 0xff, .green = 0xff, .blue = 0xff},
    .text = {.red = 0x00, .green = 0x00, .blue = 0x00},
    .fails = {.red = 0xd4, .green = 0xd4, .blue = 0xd4},
    .meets = {.red = 0x8c, .green = 0xcf, .blue = 0xf0},
    .fails_border = {.red = 0x6b, .green = 0x6b, .blue = 0x6b},
    .meets_border = {.red = 0x00, .green = 0x5a, .blue = 0x8c},
};

constexpr auto dark_picking_palette = PickingPalette{
    .view_background = {.red = 0x1e, .green = 0x1e, .blue = 0x1e},
    .text = {.red = 0xff, .green = 0xff, .blue = 0xff},
    .fails = {.red = 0x55, .green = 0x55, .blue = 0x55},
    .meets = {.red = 0x00, .green = 0x5f, .blue = 0x8f},
    .fails_border = {.red = 0xa6, .green = 0xa6, .blue = 0xa6},
    .meets_border = {.red = 0x7f, .green = 0xd0, .blue = 0xf5},
};

}  // namespace primerpicker::gui
