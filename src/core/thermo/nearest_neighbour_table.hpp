#pragma once

#include "core/sequence/nucleotide.hpp"

#include <array>
#include <cstddef>
#include <iterator>

namespace primerpicker::thermo {

// A table of nearest-neighbour values, indexed by [5' nucleotide][3'
// nucleotide] of a pair on one strand, in Nucleotide order (A, C, G, T)
template <typename Value>
using NearestNeighbourTable = std::array<std::array<Value, 4>, 4>;

template <typename Value>
[[nodiscard]] constexpr auto lookup(NearestNeighbourTable<Value> const& table, sequence::Nucleotide const first,
                                    sequence::Nucleotide const second) noexcept -> Value const& {
    // the enum has four values, so both indexes are within the table
    auto const row = std::next(table.cbegin(), static_cast<std::ptrdiff_t>(first));
    return *std::next(row->cbegin(), static_cast<std::ptrdiff_t>(second));
}

}  // namespace primerpicker::thermo
