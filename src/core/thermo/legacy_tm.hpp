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

}  // namespace primerpicker::thermo
