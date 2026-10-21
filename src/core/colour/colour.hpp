#pragma once

#include <cstdint>

namespace primerpicker::colour {

// A colour of the sRGB colour space, 8 bits per component
struct Srgb {
    std::uint8_t red;
    std::uint8_t green;
    std::uint8_t blue;
};

// The same, components from 0 to 1 (gamma-encoded, as Srgb / 255)
struct SrgbComponents {
    double red;
    double green;
    double blue;
};

[[nodiscard]] constexpr auto components(Srgb const colour) noexcept -> SrgbComponents {
    constexpr auto maximum = 255.0;
    return {.red = colour.red / maximum, .green = colour.green / maximum, .blue = colour.blue / maximum};
}

// Relative luminance, from 0 (black) to 1 (white): WCAG 2.2, definition
// of "relative luminance" (sRGB)
[[nodiscard]] auto relative_luminance(SrgbComponents colour) noexcept -> double;
[[nodiscard]] auto relative_luminance(Srgb colour) noexcept -> double;

// Contrast ratio of two colours, from 1 to 21 (WCAG 2.2, definition of
// "contrast ratio"); the order of the colours does not matter. WCAG
// AA asks for 4.5 between normal text and its background.
[[nodiscard]] auto contrast_ratio(Srgb first, Srgb second) noexcept -> double;

enum struct ColourVision : std::uint8_t { protanopia, deuteranopia, tritanopia };

// How a colour looks with a colour-vision deficiency: the simulation of
// Machado, Oliveira and Fernandes (2009, IEEE Transactions on
// Visualization and Computer Graphics 15(6):1291-1298), severity 1,
// applied to linear sRGB; components clipped to [0, 1]
[[nodiscard]] auto simulate(Srgb colour, ColourVision vision) noexcept -> SrgbComponents;

// The CIE 1976 colour difference (Delta E*ab) of two colours, in
// CIELAB with the D65 white point; about 2.3 is a just noticeable
// difference
[[nodiscard]] auto colour_difference_cie76(SrgbComponents first, SrgbComponents second) noexcept -> double;
[[nodiscard]] auto colour_difference_cie76(Srgb first, Srgb second) noexcept -> double;

}  // namespace primerpicker::colour
