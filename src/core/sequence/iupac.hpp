#pragma once

#include "core/sequence/nucleotide.hpp"

#include <cstdint>
#include <initializer_list>
#include <numeric>
#include <optional>

namespace primerpicker::sequence {

// A set of nucleotides, as an IUPAC letter stands for one (R: A or G)
class NucleotideSet {
public:
    constexpr NucleotideSet() noexcept = default;

    constexpr NucleotideSet(std::initializer_list<Nucleotide> const nucleotides) noexcept
        : bits_{std::accumulate(nucleotides.begin(), nucleotides.end(), std::uint8_t{0},
                                [](std::uint8_t const bits, Nucleotide const nucleotide) {
                                    return static_cast<std::uint8_t>(bits | bit(nucleotide));
                                })} {}

    [[nodiscard]] constexpr auto contains(Nucleotide const nucleotide) const noexcept -> bool {
        return (bits_ & bit(nucleotide)) != 0U;
    }

    // the complements of the nucleotides (R, A or G: Y, T or C)
    [[nodiscard]] constexpr auto complement() const noexcept -> NucleotideSet {
        auto result = NucleotideSet{};
        for (auto const nucleotide :
             {Nucleotide::adenine, Nucleotide::cytosine, Nucleotide::guanine, Nucleotide::thymine}) {
            if (contains(nucleotide)) {
                result.bits_ |= bit(sequence::complement(nucleotide));
            }
        }
        return result;
    }

    friend constexpr auto operator==(NucleotideSet const&, NucleotideSet const&) noexcept -> bool = default;

private:
    [[nodiscard]] static constexpr auto bit(Nucleotide const nucleotide) noexcept -> std::uint8_t {
        return static_cast<std::uint8_t>(1U << static_cast<unsigned>(nucleotide));
    }

    std::uint8_t bits_ = 0;
};

// The nucleotides an IUPAC letter stands for (NC-IUB 1985, Eur J
// Biochem 150:1-5, Table 1), in either case; other letters have no value
[[nodiscard]] constexpr auto from_iupac(char const letter) noexcept -> std::optional<NucleotideSet> {
    using enum Nucleotide;
    constexpr auto lower_to_upper = 'a' - 'A';
    auto const upper = letter >= 'a' and letter <= 'z' ? static_cast<char>(letter - lower_to_upper) : letter;
    switch (upper) {
        case 'A':
            return NucleotideSet{adenine};
        case 'C':
            return NucleotideSet{cytosine};
        case 'G':
            return NucleotideSet{guanine};
        case 'T':
            return NucleotideSet{thymine};
        case 'R':
            return NucleotideSet{adenine, guanine};
        case 'Y':
            return NucleotideSet{cytosine, thymine};
        case 'S':
            return NucleotideSet{cytosine, guanine};
        case 'W':
            return NucleotideSet{adenine, thymine};
        case 'K':
            return NucleotideSet{guanine, thymine};
        case 'M':
            return NucleotideSet{adenine, cytosine};
        case 'B':
            return NucleotideSet{cytosine, guanine, thymine};
        case 'D':
            return NucleotideSet{adenine, guanine, thymine};
        case 'H':
            return NucleotideSet{adenine, cytosine, thymine};
        case 'V':
            return NucleotideSet{adenine, cytosine, guanine};
        case 'N':
            return NucleotideSet{adenine, cytosine, guanine, thymine};
        default:
            return std::nullopt;
    }
}

}  // namespace primerpicker::sequence
