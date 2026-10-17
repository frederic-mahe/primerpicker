#pragma once

#include "core/sequence/nucleotide.hpp"
#include "core/thermo/units.hpp"

#include <optional>
#include <span>

namespace primerpicker::thermo {

// Primer Tm of Primer Designer 2.0, reverse-engineered from recorded
// values (tests/legacy/README.md, "Primer Tm"):
//   Tm = 1000 dH / (dS + b) - 273
// with Breslauer 1986 sums (breslauer_1986_sums). The 5' extension of
// a primer is not part of its Tm. Requires at least two nucleotides.
[[nodiscard]] auto legacy_primer_tm(std::span<sequence::Nucleotide const> primer) noexcept
    -> Celsius;

// The value Primer Designer shows: truncated, and none ("?" on screen)
// when not positive
[[nodiscard]] auto legacy_displayed_tm(Celsius temperature) noexcept -> std::optional<int>;

// Product Tm of Primer Designer 2.0 (tests/legacy/README.md, "Product
// Tm and annealing temperature"):
//   Tm = 59.9 + 41 x (GC fraction) - 675 / length
// that is 81.5 + 16.6 log10(0.05) + 0.41 %GC - 675 / length. Requires
// a non-empty product.
[[nodiscard]] auto legacy_product_tm(std::span<sequence::Nucleotide const> product) noexcept
    -> Celsius;

// The displayed (truncated) Tm of both primers and of their product
struct LegacyPairTms {
    int primer_a_celsius;
    int primer_b_celsius;
    int product_celsius;
};

// Annealing temperature of Primer Designer 2.0, after Rychlik et al.
// 1990 (NAR 18:6409), from the displayed Tm values:
//   Ta = 0.3 x (lower primer Tm - 21.6) + 0.7 x product Tm - 14.9
[[nodiscard]] auto legacy_annealing_tm(LegacyPairTms const& tms) noexcept -> Celsius;

}  // namespace primerpicker::thermo
