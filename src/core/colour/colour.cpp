#include "core/colour/colour.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

// C++26 refactoring: make these functions constexpr (std::pow, std::cbrt)

namespace primerpicker::colour {

namespace {

using Matrix = std::array<std::array<double, 3>, 3>;

struct Vector {
    double first;
    double second;
    double third;
};

constexpr auto multiply(Matrix const& matrix, Vector const& vector) noexcept -> Vector {
    auto const row = [&vector](std::array<double, 3> const& coefficients) {
        return (coefficients[0] * vector.first) + (coefficients[1] * vector.second) + (coefficients[2] * vector.third);
    };
    return {.first = row(matrix[0]), .second = row(matrix[1]), .third = row(matrix[2])};
}

// sRGB transfer functions (IEC 61966-2-1)
auto to_linear(double const component) noexcept -> double {
    constexpr auto threshold = 0.04045;
    constexpr auto slope = 12.92;
    constexpr auto offset = 0.055;
    constexpr auto exponent = 2.4;
    if (component <= threshold) {
        return component / slope;
    }
    return std::pow((component + offset) / (1.0 + offset), exponent);
}

auto from_linear(double const component) noexcept -> double {
    constexpr auto threshold = 0.0031308;
    constexpr auto slope = 12.92;
    constexpr auto offset = 0.055;
    constexpr auto exponent = 1.0 / 2.4;
    if (component <= threshold) {
        return component * slope;
    }
    return ((1.0 + offset) * std::pow(component, exponent)) - offset;
}

auto linear(SrgbComponents const colour) noexcept -> Vector {
    return {.first = to_linear(colour.red), .second = to_linear(colour.green), .third = to_linear(colour.blue)};
}

// Machado, Oliveira and Fernandes 2009: the matrices published with the
// paper, severity 1.0; rows give the simulated linear R, G and B
constexpr auto protanopia_matrix = Matrix{{{0.152286, 1.052583, -0.204868},
                                           {0.114503, 0.786281, 0.099216},
                                           {-0.003882, -0.048116, 1.051998}}};
constexpr auto deuteranopia_matrix = Matrix{{{0.367322, 0.860646, -0.227968},
                                             {0.280085, 0.672501, 0.047413},
                                             {-0.011820, 0.042940, 0.968881}}};
constexpr auto tritanopia_matrix = Matrix{{{1.255528, -0.076749, -0.178779},
                                           {-0.078411, 0.930809, 0.147602},
                                           {0.004733, 0.691367, 0.303900}}};

constexpr auto vision_matrix(ColourVision const vision) noexcept -> Matrix const& {
    switch (vision) {
        case ColourVision::protanopia:
            return protanopia_matrix;
        case ColourVision::deuteranopia:
            return deuteranopia_matrix;
        case ColourVision::tritanopia:
            return tritanopia_matrix;
        default:
            return protanopia_matrix;  // unreachable: all enumerators are handled
    }
}

// CIE XYZ (Y of white 100) to linear sRGB: the matrix of IEC
// 61966-2-1:1999, whose inverse gives XYZ (as colorspacious 1.1.2 does)
constexpr auto srgb_matrix = Matrix{{{3.2406, -1.5372, -0.4986},  //
                                     {-0.9689, 1.8758, 0.0415},
                                     {0.0557, -0.2040, 1.0570}}};

constexpr auto inverse(Matrix const& matrix) noexcept -> Matrix {
    auto const& first = matrix[0];
    auto const cofactor = [&matrix](std::size_t const row, std::size_t const column) {
        auto const next_row = (row + 1) % 3;
        auto const last_row = (row + 2) % 3;
        auto const next_column = (column + 1) % 3;
        auto const last_column = (column + 2) % 3;
        return (matrix.at(next_row).at(next_column) * matrix.at(last_row).at(last_column)) -
               (matrix.at(next_row).at(last_column) * matrix.at(last_row).at(next_column));
    };
    auto const determinant =
        (first[0] * cofactor(0, 0)) + (first[1] * cofactor(0, 1)) + (first[2] * cofactor(0, 2));
    auto result = Matrix{};
    for (auto row = std::size_t{0}; row < 3; ++row) {
        for (auto column = std::size_t{0}; column < 3; ++column) {
            result.at(row).at(column) = cofactor(column, row) / determinant;  // the adjugate, transposed
        }
    }
    return result;
}

constexpr auto xyz_matrix = inverse(srgb_matrix);  // linear sRGB to XYZ, Y of white 1
// the D65 white point, Y = 1 (colorspacious 1.1.2's: X 95.047 and
// Z 108.883 for Y 100)
constexpr auto white_point = Vector{.first = 0.95047, .second = 1.0, .third = 1.08883};

struct Lab {
    double lightness;
    double a_axis;
    double b_axis;
};

// CIELAB (CIE 1976): f(t) of the XYZ ratios to the white point
auto lab_function(double const ratio) noexcept -> double {
    constexpr auto delta = 6.0 / 29.0;
    if (ratio > delta * delta * delta) {
        return std::cbrt(ratio);
    }
    constexpr auto four_29ths = 4.0 / 29.0;
    return (ratio / (3.0 * delta * delta)) + four_29ths;
}

auto lab(SrgbComponents const colour) noexcept -> Lab {
    auto const xyz = multiply(xyz_matrix, linear(colour));
    auto const fx = lab_function(xyz.first / white_point.first);
    auto const fy = lab_function(xyz.second / white_point.second);
    auto const fz = lab_function(xyz.third / white_point.third);
    constexpr auto lightness_scale = 116.0;
    constexpr auto lightness_offset = 16.0;
    constexpr auto a_scale = 500.0;
    constexpr auto b_scale = 200.0;
    return {.lightness = (lightness_scale * fy) - lightness_offset,
            .a_axis = a_scale * (fx - fy),
            .b_axis = b_scale * (fy - fz)};
}

}  // namespace

auto relative_luminance(SrgbComponents const colour) noexcept -> double {
    // WCAG 2.2 coefficients (those of ITU-R BT.709)
    constexpr auto red_weight = 0.2126;
    constexpr auto green_weight = 0.7152;
    constexpr auto blue_weight = 0.0722;
    // WCAG 2.2 writes 0.03928 for the threshold; with 8-bit components
    // both thresholds give the same values
    auto const [red, green, blue] = linear(colour);
    return (red_weight * red) + (green_weight * green) + (blue_weight * blue);
}

auto relative_luminance(Srgb const colour) noexcept -> double {
    return relative_luminance(components(colour));
}

auto contrast_ratio(Srgb const first, Srgb const second) noexcept -> double {
    constexpr auto flare = 0.05;
    auto const one = relative_luminance(first);
    auto const other = relative_luminance(second);
    return (std::max(one, other) + flare) / (std::min(one, other) + flare);
}

auto simulate(Srgb const colour, ColourVision const vision) noexcept -> SrgbComponents {
    auto const simulated = multiply(vision_matrix(vision), linear(components(colour)));
    auto const encoded = [](double const component) { return from_linear(std::clamp(component, 0.0, 1.0)); };
    return {.red = encoded(simulated.first), .green = encoded(simulated.second), .blue = encoded(simulated.third)};
}

auto colour_difference_cie76(SrgbComponents const first, SrgbComponents const second) noexcept -> double {
    auto const one = lab(first);
    auto const other = lab(second);
    return std::hypot(one.lightness - other.lightness, one.a_axis - other.a_axis, one.b_axis - other.b_axis);
}

auto colour_difference_cie76(Srgb const first, Srgb const second) noexcept -> double {
    return colour_difference_cie76(components(first), components(second));
}

}  // namespace primerpicker::colour
