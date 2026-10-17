#pragma once

#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

namespace primerpicker::sequence {

enum struct Nucleotide : std::uint8_t { adenine, cytosine, guanine, thymine };

// A, C, G or T, in either case; anything else has no value
[[nodiscard]] auto to_nucleotide(char letter) noexcept -> std::optional<Nucleotide>;

// no value if any letter is not a nucleotide
[[nodiscard]] auto to_nucleotides(std::string_view letters)
    -> std::optional<std::vector<Nucleotide>>;

}  // namespace primerpicker::sequence
