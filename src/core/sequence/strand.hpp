#pragma once

#include "core/sequence/nucleotide.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <optional>
#include <span>
#include <vector>

namespace primerpicker::sequence {

enum struct Strand : std::uint8_t { top, complement };

// Where a primer lies on a molecule, as Primer Designer gives it: the
// 1-based top-strand coordinate of the primer's 5' end; on the
// complement strand the primer runs towards lower coordinates
// (position 400(C), 20 bases: the complement of 381 to 400)
struct PrimerLocation {
    std::size_t position;
    std::size_t length;
    Strand strand;
};

// The last position of a search window that starts at target and
// extends offset bases beyond it, within a molecule of molecule_size
// bases (no primer lies beyond it); never overflows
[[nodiscard]] constexpr auto window_end(std::size_t const target, std::size_t const offset,
                                        std::size_t const molecule_size) noexcept -> std::size_t {
    if (target >= molecule_size or offset >= molecule_size - target) {
        return molecule_size;
    }
    return target + offset;
}

// The primer, 5' to 3'; none if it does not fit on the molecule. Not
// noexcept: the vector allocates.
[[nodiscard]] constexpr auto primer_at(std::span<Nucleotide const> const molecule, PrimerLocation const& location)
    -> std::optional<std::vector<Nucleotide>> {
    auto const& [position, length, strand] = location;
    if (position == 0 or length == 0) {
        return std::nullopt;
    }
    if (strand == Strand::top) {
        if (position > molecule.size() or length > molecule.size() - (position - 1)) {  // no overflow
            return std::nullopt;
        }
        auto const stretch = molecule.subspan(position - 1, length);
        return std::vector<Nucleotide>(stretch.begin(), stretch.end());
    }
    if (position > molecule.size() or position < length) {
        return std::nullopt;
    }
    auto const stretch = molecule.subspan(position - length, length);
    auto primer = std::vector<Nucleotide>{};
    primer.reserve(length);
    std::ranges::transform(stretch.rbegin(), stretch.rend(), std::back_inserter(primer), complement);
    return primer;
}

}  // namespace primerpicker::sequence
