#include "core/colour/colour.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace {

using Catch::Matchers::WithinAbs;
using primerpicker::colour::colour_difference_cie76;
using primerpicker::colour::ColourVision;
using primerpicker::colour::contrast_ratio;
using primerpicker::colour::SrgbComponents;
using primerpicker::colour::relative_luminance;
using primerpicker::colour::simulate;
using primerpicker::colour::Srgb;

auto parse_hex_colour(std::string_view const text) -> Srgb {
    constexpr auto hex_length = std::size_t{7};  // #rrggbb
    constexpr auto base = 16;
    if (text.size() != hex_length or not text.starts_with('#')) {
        throw std::runtime_error{"not a #rrggbb colour: " + std::string{text}};
    }
    auto components = std::array<std::uint8_t, 3>{};
    for (auto index = std::size_t{0}; index < components.size(); ++index) {
        auto const digits = text.substr(1 + (2 * index), 2);
        auto const [end, error] = std::from_chars(digits.begin(), digits.end(), components.at(index), base);
        if (error != std::errc{} or end != digits.end()) {
            throw std::runtime_error{"not a #rrggbb colour: " + std::string{text}};
        }
    }
    return {.red = components[0], .green = components[1], .blue = components[2]};
}

auto parse_components(std::string const& text) -> SrgbComponents {
    auto stream = std::istringstream{text};
    auto components = SrgbComponents{};
    auto comma = ',';
    if (not(stream >> components.red >> comma >> components.green >> comma >> components.blue)) {
        throw std::runtime_error{"not three components: " + text};
    }
    return components;
}

struct Reference {
    Srgb first;
    Srgb second;
    double luminance;
    double contrast;
    std::array<SrgbComponents, 3> simulations;  // protanopia, deuteranopia, tritanopia
    double difference;
};

// tests/data/colorspacious_colours.tsv (see colorspacious_colours.py)
auto references() -> std::vector<Reference> {
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_TEST_DATA_DIR} / "colorspacious_colours.tsv"};
    auto found = std::vector<Reference>{};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with("# ")) {  // the header; colours start with '#' too
            continue;
        }
        auto stream = std::istringstream{line};
        auto first = std::string{};
        auto second = std::string{};
        auto luminance = 0.0;
        auto contrast = 0.0;
        auto protanopia = std::string{};
        auto deuteranopia = std::string{};
        auto tritanopia = std::string{};
        auto difference = 0.0;
        if (not(stream >> first >> second >> luminance >> contrast >> protanopia >> deuteranopia >> tritanopia >>
                difference)) {
            throw std::runtime_error{"cannot read line: " + line};
        }
        found.push_back({.first = parse_hex_colour(first),
                         .second = parse_hex_colour(second),
                         .luminance = luminance,
                         .contrast = contrast,
                         .simulations = {parse_components(protanopia), parse_components(deuteranopia),
                                         parse_components(tritanopia)},
                         .difference = difference});
    }
    return found;
}

constexpr auto black = Srgb{.red = 0, .green = 0, .blue = 0};
constexpr auto white = Srgb{.red = 255, .green = 255, .blue = 255};

}  // namespace

// Reference: WCAG 2.2, definition of "contrast ratio" (ratios range
// from 1:1 to 21:1, black on white)
TEST_CASE("black on white has the largest contrast ratio, 21:1", "[colour]") {
    CHECK_THAT(contrast_ratio(black, white), WithinAbs(21.0, 1e-9));
    CHECK_THAT(contrast_ratio(white, black), WithinAbs(21.0, 1e-9));
    CHECK_THAT(contrast_ratio(white, white), WithinAbs(1.0, 1e-9));
}

// Reference: WebAIM, Contrast Checker: #767676 on white is 4.54:1, the
// lightest grey that passes WCAG AA for normal text (4.5:1); #777777
// is 4.48:1
TEST_CASE("the lightest grey that passes WCAG AA on white", "[colour]") {
    constexpr auto passes = Srgb{.red = 0x76, .green = 0x76, .blue = 0x76};
    constexpr auto fails = Srgb{.red = 0x77, .green = 0x77, .blue = 0x77};
    constexpr auto tolerance = 0.005;  // WebAIM prints two decimals
    CHECK_THAT(contrast_ratio(passes, white), WithinAbs(4.54, tolerance));
    CHECK_THAT(contrast_ratio(fails, white), WithinAbs(4.48, tolerance));
}

// Reference values: colorspacious 1.1.2 (Machado et al. 2009 matrices,
// severity 100; deltaE in CIELab) and the WCAG 2.2 formulas, written
// by tests/data/colorspacious_colours.py
TEST_CASE("luminance, contrast, colour-vision simulations and CIE76 differences", "[colour][colorspacious]") {
    auto const found = references();
    REQUIRE(found.size() == 11);
    constexpr auto tolerance = 1e-5;  // the table has 6 decimals
    constexpr auto kinds = std::array{ColourVision::protanopia, ColourVision::deuteranopia, ColourVision::tritanopia};
    for (auto const& reference : found) {
        INFO(static_cast<int>(reference.first.red) << ',' << static_cast<int>(reference.first.green) << ','
                                                   << static_cast<int>(reference.first.blue));
        CHECK_THAT(relative_luminance(reference.first), WithinAbs(reference.luminance, tolerance));
        CHECK_THAT(contrast_ratio(reference.first, reference.second), WithinAbs(reference.contrast, tolerance));
        for (auto index = std::size_t{0}; index < kinds.size(); ++index) {
            auto const simulated = simulate(reference.first, kinds.at(index));
            auto const& expected = reference.simulations.at(index);
            CHECK_THAT(simulated.red, WithinAbs(expected.red, tolerance));
            CHECK_THAT(simulated.green, WithinAbs(expected.green, tolerance));
            CHECK_THAT(simulated.blue, WithinAbs(expected.blue, tolerance));
        }
        CHECK_THAT(colour_difference_cie76(reference.first, reference.second),
                   WithinAbs(reference.difference, tolerance));
    }
}
