#pragma once

#include "core/sequence/nucleotide.hpp"

#include <cstddef>
#include <string_view>

namespace primerpicker::sequence {

namespace detail {

// The standard genetic code: NCBI, "The Genetic Codes", translation
// table 1. The amino acids (one-letter codes) of the 64 codons, the
// three bases running through T, C, A, G in turn (first base slowest);
// '-' is a stop codon, as Primer Designer 2.0 shows it
constexpr auto standard_code = std::string_view{"FFLLSSSSYY--CC-WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"};
constexpr auto bases_per_codon_position = std::size_t{4};

// the rank of a base in the T, C, A, G order of the table
constexpr auto tcag_rank(Nucleotide const nucleotide) noexcept -> std::size_t {
    switch (nucleotide) {
        case Nucleotide::thymine:
            return 0;
        case Nucleotide::cytosine:
            return 1;
        case Nucleotide::adenine:
            return 2;
        case Nucleotide::guanine:
            return 3;
        default:
            return 0;  // unreachable: all enumerators are handled
    }
}

}  // namespace detail

// The amino acid (one-letter code) of a codon, '-' for a stop codon
[[nodiscard]] constexpr auto translate(Nucleotide const first, Nucleotide const second,
                                       Nucleotide const third) noexcept -> char {
    constexpr auto base = detail::bases_per_codon_position;
    auto const index = (detail::tcag_rank(first) * base * base) + (detail::tcag_rank(second) * base) +
                       detail::tcag_rank(third);
    return detail::standard_code[index];
}

}  // namespace primerpicker::sequence
