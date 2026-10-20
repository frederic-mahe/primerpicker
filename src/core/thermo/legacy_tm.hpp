#pragma once

#include "core/sequence/nucleotide.hpp"
#include "core/thermo/breslauer_1986.hpp"
#include "core/thermo/units.hpp"

#include <algorithm>
#include <cassert>
#include <optional>
#include <span>

namespace primerpicker::thermo {

namespace detail {

// Fitted on 186 primers recorded from PRIMER.exe: every displayed
// value matches for b in [60.6983, 60.7000] cal/K/mol
// (tests/legacy/tools/fit_tm.py). The origin of b is not known.
// Primers recorded on 2026-10-20 narrow it to (60.699732, 60.7]
// (tests/legacy/README.md, "Primer Tm: the constant b"); 60.7 is
// probably the original value, but the middle of the interval
// (60.69987) is used, so that a Tm computed exactly on an integer is
// not truncated below it.
constexpr auto legacy_entropy_constant = CaloriesPerKelvinPerMole{60.69987};
// Primer3's oligotm gives the same Tm in Breslauer mode with a DNA
// concentration C of 49.6 pM: b = 10.8 + R ln(4 / C), 10.8 being
// Primer3's initiation entropy (tests/data/primer3_legacy_tm.sh)

// Primer Designer subtracts 273, not 273.15: 273.15 contradicts
// several recorded values
constexpr auto legacy_kelvin_offset = 273.0;

constexpr auto cal_per_kcal = 1000.0;

// Product Tm and annealing temperature: these doubles are stored next
// to each other in PRIMER.exe (file offset 0x429d3), and match 28
// recorded pairs (tests/legacy/README.md)
constexpr auto product_tm_base = 59.9;  // 81.5 + 16.6 log10(0.05)
constexpr auto product_tm_gc_factor = 41.0;  // 0.41 x %GC
constexpr auto product_tm_length_factor = 675.0;
constexpr auto annealing_primer_offset = 21.6;  // -16.6 log10(0.05)
constexpr auto annealing_primer_weight = 0.3;
constexpr auto annealing_product_weight = 0.7;
constexpr auto annealing_offset = 14.9;

}  // namespace detail

// Primer Tm of Primer Designer 2.0, reverse-engineered from recorded
// values (tests/legacy/README.md, "Primer Tm"):
//   Tm = 1000 dH / (dS + b) - 273
// with Breslauer 1986 sums (breslauer_1986_sums). The 5' extension of
// a primer is not part of its Tm. Requires at least two nucleotides.
[[nodiscard]] constexpr auto legacy_primer_tm(std::span<sequence::Nucleotide const> const primer) noexcept
    -> Celsius {
    assert(primer.size() >= 2);
    auto const sums = breslauer_1986_sums(primer);
    auto const entropy =
        sums.entropy.cal_per_kelvin_per_mol + detail::legacy_entropy_constant.cal_per_kelvin_per_mol;
    return Celsius{(detail::cal_per_kcal * sums.enthalpy.kcal_per_mol / entropy) - detail::legacy_kelvin_offset};
}

// The value Primer Designer shows: truncated, and none ("?" on screen)
// when not positive
// C++23 refactoring: make constexpr (std::floor, std::isfinite)
[[nodiscard]] auto legacy_displayed_tm(Celsius temperature) noexcept -> std::optional<int>;

// Product Tm of Primer Designer 2.0 (tests/legacy/README.md, "Product
// Tm and annealing temperature"):
//   Tm = 59.9 + 41 x (GC fraction) - 675 / length
// that is 81.5 + 16.6 log10(0.05) + 0.41 %GC - 675 / length. Requires
// a non-empty product.
[[nodiscard]] constexpr auto legacy_product_tm(std::span<sequence::Nucleotide const> const product) noexcept
    -> Celsius {
    assert(not product.empty());
    auto const is_gc = [](sequence::Nucleotide const nucleotide) {
        return nucleotide == sequence::Nucleotide::cytosine or nucleotide == sequence::Nucleotide::guanine;
    };
    auto const length = static_cast<double>(product.size());
    auto const gc_fraction = static_cast<double>(std::ranges::count_if(product, is_gc)) / length;
    return Celsius{detail::product_tm_base + (detail::product_tm_gc_factor * gc_fraction) -
                   (detail::product_tm_length_factor / length)};
}

// The displayed (truncated) Tm of both primers and of their product
struct LegacyPairTms {
    int primer_a_celsius;
    int primer_b_celsius;
    int product_celsius;
};

// Annealing temperature of Primer Designer 2.0, after Rychlik et al.
// 1990 (NAR 18:6409), from the displayed Tm values:
//   Ta = 0.3 x (lower primer Tm - 21.6) + 0.7 x product Tm - 14.9
[[nodiscard]] constexpr auto legacy_annealing_tm(LegacyPairTms const& tms) noexcept -> Celsius {
    auto const lower_primer_tm = std::min(tms.primer_a_celsius, tms.primer_b_celsius);
    return Celsius{(detail::annealing_primer_weight * (lower_primer_tm - detail::annealing_primer_offset)) +
                   (detail::annealing_product_weight * tms.product_celsius) - detail::annealing_offset};
}

}  // namespace primerpicker::thermo
