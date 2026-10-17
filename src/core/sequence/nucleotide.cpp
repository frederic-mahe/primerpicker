#include "core/sequence/nucleotide.hpp"

#include <optional>
#include <string_view>
#include <vector>

namespace primerpicker::sequence {

auto to_nucleotide(char const letter) noexcept -> std::optional<Nucleotide> {
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

// not noexcept: the vector allocates
auto to_nucleotides(std::string_view const letters)
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
