#pragma once

#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

namespace primerpicker::sequence {

enum struct Nucleotide : std::uint8_t { adenine, cytosine, guanine, thymine };

// Watson-Crick complement
[[nodiscard]] constexpr auto complement(Nucleotide const nucleotide) noexcept -> Nucleotide {
    switch (nucleotide) {
        case Nucleotide::adenine:
            return Nucleotide::thymine;
        case Nucleotide::cytosine:
            return Nucleotide::guanine;
        case Nucleotide::guanine:
            return Nucleotide::cytosine;
        case Nucleotide::thymine:
            return Nucleotide::adenine;
        default:
            return nucleotide;  // unreachable: all enumerators are handled
    }
}

// The upper-case letter of a nucleotide
[[nodiscard]] constexpr auto to_letter(Nucleotide const nucleotide) noexcept -> char {
    switch (nucleotide) {
        case Nucleotide::adenine:
            return 'A';
        case Nucleotide::cytosine:
            return 'C';
        case Nucleotide::guanine:
            return 'G';
        case Nucleotide::thymine:
            return 'T';
        default:
            return '?';  // unreachable: all enumerators are handled
    }
}

// A, C, G or T, in either case; anything else has no value
[[nodiscard]] constexpr auto to_nucleotide(char const letter) noexcept -> std::optional<Nucleotide> {
    switch (letter) {
        case 'A':
        case 'a':
            return Nucleotide::adenine;
        case 'C':
        case 'c':
            return Nucleotide::cytosine;
        case 'G':
        case 'g':
            return Nucleotide::guanine;
        case 'T':
        case 't':
            return Nucleotide::thymine;
        default:
            return std::nullopt;
    }
}

// no value if any letter is not a nucleotide; not noexcept: the vector
// allocates
[[nodiscard]] constexpr auto to_nucleotides(std::string_view const letters)
    -> std::optional<std::vector<Nucleotide>> {
    auto nucleotides = std::vector<Nucleotide>{};
    nucleotides.reserve(letters.size());
    for (auto const letter : letters) {
        auto const nucleotide = to_nucleotide(letter);
        if (not nucleotide.has_value()) {
            return std::nullopt;
        }
        nucleotides.push_back(*nucleotide);
    }
    return nucleotides;
}

}  // namespace primerpicker::sequence
