#pragma once

#include "core/sequence/nucleotide.hpp"

#include <cstddef>

namespace primerpicker::sequence {

// The Direct Entry field of Primer Designer 2.0 (Primer/Create ->
// Direct Entry; tests/legacy/README.md, "Primer Tm"): it holds 50
// bases, further keys are ignored
constexpr auto direct_entry_max_length = std::size_t{50};

// Whether the Direct Entry field inserts a typed letter: A, C, G and T
// only (PRIMER.exe ignores N, R and U). Lower-case letters are not
// recorded; they are inserted, as the command line accepts them.
[[nodiscard]] constexpr auto direct_entry_inserts(char const letter) noexcept -> bool {
    return to_nucleotide(letter).has_value();
}

}  // namespace primerpicker::sequence
