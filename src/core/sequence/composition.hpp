#pragma once

#include "core/sequence/nucleotide.hpp"

#include <algorithm>
#include <cstddef>
#include <span>

namespace primerpicker::sequence {

// The number of each nucleotide in a sequence ("Composition: A: C: G:
// T:" of Primer Designer's screens and reports)
struct Composition {
    std::size_t adenine = 0;
    std::size_t cytosine = 0;
    std::size_t guanine = 0;
    std::size_t thymine = 0;

    friend constexpr auto operator==(Composition const&, Composition const&) noexcept -> bool = default;
};

[[nodiscard]] constexpr auto composition(std::span<Nucleotide const> const sequence) noexcept -> Composition {
    auto counts = Composition{};
    std::ranges::for_each(sequence, [&counts](Nucleotide const nucleotide) {
        switch (nucleotide) {
            case Nucleotide::adenine:
                ++counts.adenine;
                break;
            case Nucleotide::cytosine:
                ++counts.cytosine;
                break;
            case Nucleotide::guanine:
                ++counts.guanine;
                break;
            case Nucleotide::thymine:
                ++counts.thymine;
                break;
            default:
                break;  // unreachable: all enumerators are handled
        }
    });
    return counts;
}

}  // namespace primerpicker::sequence
