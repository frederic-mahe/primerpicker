#include "core/colour/colour.hpp"
#include "gui/picking_colours.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>

namespace {

using primerpicker::colour::colour_difference_cie76;
using primerpicker::colour::ColourVision;
using primerpicker::colour::components;
using primerpicker::colour::contrast_ratio;
using primerpicker::colour::simulate;
using primerpicker::gui::dark_picking_palette;
using primerpicker::gui::light_picking_palette;
using primerpicker::gui::PickingPalette;

constexpr auto palettes = std::array{light_picking_palette, dark_picking_palette};

// WCAG 2.2, success criterion 1.4.3 (AA): 4.5:1 for normal text
constexpr auto text_contrast = 4.5;
// WCAG 2.2, success criterion 1.4.11: 3:1 for user interface
// components and graphical objects (the window's border)
constexpr auto border_contrast = 3.0;
// The two window backgrounds must stay clearly apart, with every colour
// vision: a CIE76 difference of 20 or more (about 9 just noticeable
// differences; the project's own threshold, proposal §6.3e)
constexpr auto minimum_difference = 20.0;

}  // namespace

TEST_CASE("the bases stay readable on the picking window, in both schemes (WCAG AA)", "[gui][picking][colour]") {
    for (auto const& palette : palettes) {
        CHECK(contrast_ratio(palette.text, palette.view_background) >= text_contrast);
        CHECK(contrast_ratio(palette.text, palette.fails) >= text_contrast);
        CHECK(contrast_ratio(palette.text, palette.meets) >= text_contrast);
    }
}

TEST_CASE("the window's borders stand out from the page and from the window", "[gui][picking][colour]") {
    for (auto const& palette : palettes) {
        CHECK(contrast_ratio(palette.fails_border, palette.view_background) >= border_contrast);
        CHECK(contrast_ratio(palette.meets_border, palette.view_background) >= border_contrast);
        CHECK(contrast_ratio(palette.fails_border, palette.fails) >= border_contrast);
        CHECK(contrast_ratio(palette.meets_border, palette.meets) >= border_contrast);
    }
}

TEST_CASE("'meets' and 'fails' differ with protanopia, deuteranopia and tritanopia too", "[gui][picking][colour]") {
    constexpr auto visions = std::array{ColourVision::protanopia, ColourVision::deuteranopia, ColourVision::tritanopia};
    for (auto const& palette : palettes) {
        CHECK(colour_difference_cie76(palette.fails, palette.meets) >= minimum_difference);
        for (auto const vision : visions) {
            CHECK(colour_difference_cie76(simulate(palette.fails, vision), simulate(palette.meets, vision)) >=
                  minimum_difference);
        }
    }
}

TEST_CASE("the light palette is light and the dark palette dark", "[gui][picking][colour]") {
    auto const lighter = [](auto const& one, auto const& other) {
        return primerpicker::colour::relative_luminance(components(one)) >
               primerpicker::colour::relative_luminance(components(other));
    };
    CHECK(lighter(light_picking_palette.view_background, light_picking_palette.text));
    CHECK(lighter(dark_picking_palette.text, dark_picking_palette.view_background));
}
