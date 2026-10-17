#include "core/thermo/breslauer_1986.hpp"

#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>

namespace primerpicker::thermo {

namespace {

auto pair_parameters(sequence::Nucleotide const first, sequence::Nucleotide const second) noexcept
    -> NearestNeighbourSums {
    auto const& parameters =
        breslauer_1986_table.at(static_cast<std::size_t>(first)).at(static_cast<std::size_t>(second));
    return {KilocaloriesPerMole{parameters.enthalpy_kcal_per_mol},
            CaloriesPerKelvinPerMole{parameters.entropy_cal_per_kelvin_per_mol}};
}

auto add(NearestNeighbourSums const& lhs, NearestNeighbourSums const& rhs) noexcept
    -> NearestNeighbourSums {
    return {KilocaloriesPerMole{lhs.enthalpy.kcal_per_mol + rhs.enthalpy.kcal_per_mol},
            CaloriesPerKelvinPerMole{lhs.entropy.cal_per_kelvin_per_mol +
                                     rhs.entropy.cal_per_kelvin_per_mol}};
}

}  // namespace

auto breslauer_1986_sums(std::span<sequence::Nucleotide const> const strand) noexcept
    -> NearestNeighbourSums {
    auto const zero = NearestNeighbourSums{KilocaloriesPerMole{0.0}, CaloriesPerKelvinPerMole{0.0}};
    if (strand.size() < 2) {
        return zero;
    }
    // C++23 refactoring: replace with std::views::pairwise
    return std::transform_reduce(strand.begin(), std::prev(strand.end()), std::next(strand.begin()),
                                 zero, add, pair_parameters);
}

}  // namespace primerpicker::thermo
