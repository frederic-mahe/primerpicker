#include "core/thermo/legacy_tm.hpp"

#include "core/thermo/breslauer_1986.hpp"

#include <cassert>
#include <cmath>

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

}  // namespace

auto legacy_primer_tm(std::span<sequence::Nucleotide const> const primer) noexcept -> Celsius {
    assert(primer.size() >= 2);
    auto const sums = breslauer_1986_sums(primer);
    auto const entropy = sums.entropy.cal_per_kelvin_per_mol + entropy_constant.cal_per_kelvin_per_mol;
    return Celsius{cal_per_kcal * sums.enthalpy.kcal_per_mol / entropy - kelvin_offset};
}

auto legacy_displayed_tm(Celsius const tm) noexcept -> std::optional<int> {
    assert(std::isfinite(tm.degrees));
    auto const truncated = static_cast<int>(std::floor(tm.degrees));
    if (truncated <= 0) {
        return std::nullopt;
    }
    return truncated;
}

}  // namespace primerpicker::thermo
