#pragma once

#include "core/sequence/nucleotide.hpp"
#include "core/thermo/units.hpp"

#include <array>
#include <span>

namespace primerpicker::thermo {

struct NearestNeighbourSums {
    KilocaloriesPerMole enthalpy;
    CaloriesPerKelvinPerMole entropy;
};

struct PairParameters {
    double enthalpy_kcal_per_mol;
    double entropy_cal_per_kelvin_per_mol;
};

// Breslauer KJ, Frank R, Blöcker H, Marky LA (1986) Predicting DNA
// duplex stability from the base sequence. PNAS 83:3746-3750, Table 2.
// Indexed by [5' nucleotide][3' nucleotide] of a pair on one strand,
// in Nucleotide order (A, C, G, T); dissociation values (positive).
constexpr auto breslauer_1986_table = std::array<std::array<PairParameters, 4>, 4>{{
    // AA, AC, AG, AT
    {{{9.1, 24.0}, {6.5, 17.3}, {7.8, 20.8}, {8.6, 23.9}}},
    // CA, CC, CG, CT
    {{{5.8, 12.9}, {11.0, 26.6}, {11.9, 27.8}, {7.8, 20.8}}},
    // GA, GC, GG, GT
    {{{5.6, 13.5}, {11.1, 26.7}, {11.0, 26.6}, {6.5, 17.3}}},
    // TA, TC, TG, TT
    {{{6.0, 16.9}, {5.6, 13.5}, {5.8, 12.9}, {9.1, 24.0}}},
}};

// Sums over the nearest-neighbour pairs of a single strand; zero for
// fewer than two nucleotides. No initiation or symmetry term.
[[nodiscard]] auto breslauer_1986_sums(std::span<sequence::Nucleotide const> strand) noexcept
    -> NearestNeighbourSums;

}  // namespace primerpicker::thermo
