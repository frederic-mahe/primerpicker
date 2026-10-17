#include "core/thermo/legacy_tm.hpp"

#include "core/sequence/nucleotide.hpp"
#include "core/thermo/breslauer_1986.hpp"
#include "core/thermo/units.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <optional>
#include <span>

namespace primerpicker::thermo {

namespace {

// Fitted on 186 primers recorded from PRIMER.exe: every displayed
// value matches for b in [60.6983, 60.7000] cal/K/mol
// (tests/legacy/tools/fit_tm.py). The origin of b is not known.
constexpr auto entropy_constant = CaloriesPerKelvinPerMole{60.699};

// Primer Designer subtracts 273, not 273.15: 273.15 contradicts
// several recorded values
constexpr auto kelvin_offset = 273.0;

constexpr auto cal_per_kcal = 1000.0;

// Product Tm and annealing temperature: these doubles are stored next
// to each other in PRIMER.exe (file offset 0x429d3), and match 20
// recorded pairs (tests/legacy/README.md)
constexpr auto product_tm_base = 59.9;           // 81.5 + 16.6 log10(0.05)
constexpr auto product_tm_gc_factor = 41.0;      // 0.41 x %GC
constexpr auto product_tm_length_factor = 675.0;
constexpr auto annealing_primer_offset = 21.6;   // -16.6 log10(0.05)
constexpr auto annealing_primer_weight = 0.3;
constexpr auto annealing_product_weight = 0.7;
constexpr auto annealing_offset = 14.9;

}  // namespace

auto legacy_primer_tm(std::span<sequence::Nucleotide const> const primer) noexcept -> Celsius {
    assert(primer.size() >= 2);
    auto const sums = breslauer_1986_sums(primer);
    auto const entropy = sums.entropy.cal_per_kelvin_per_mol + entropy_constant.cal_per_kelvin_per_mol;
    return Celsius{(cal_per_kcal * sums.enthalpy.kcal_per_mol / entropy) - kelvin_offset};
}

auto legacy_displayed_tm(Celsius const temperature) noexcept -> std::optional<int> {
    assert(std::isfinite(temperature.degrees));
    auto const truncated = static_cast<int>(std::floor(temperature.degrees));
    if (truncated <= 0) {
        return std::nullopt;
    }
    return truncated;
}

auto legacy_product_tm(std::span<sequence::Nucleotide const> const product) noexcept -> Celsius {
    assert(not product.empty());
    auto const is_gc = [](sequence::Nucleotide const nucleotide) {
        return nucleotide == sequence::Nucleotide::cytosine or
               nucleotide == sequence::Nucleotide::guanine;
    };
    auto const length = static_cast<double>(product.size());
    auto const gc_fraction = static_cast<double>(std::ranges::count_if(product, is_gc)) / length;
    return Celsius{product_tm_base + (product_tm_gc_factor * gc_fraction) -
                   (product_tm_length_factor / length)};
}

auto legacy_annealing_tm(LegacyPairTms const& tms) noexcept -> Celsius {
    auto const lower_primer_tm = std::min(tms.primer_a_celsius, tms.primer_b_celsius);
    return Celsius{(annealing_primer_weight * (lower_primer_tm - annealing_primer_offset)) +
                   (annealing_product_weight * tms.product_celsius) - annealing_offset};
}

}  // namespace primerpicker::thermo
