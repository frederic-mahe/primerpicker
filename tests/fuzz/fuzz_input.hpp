#pragma once

// Helpers shared by the libFuzzer harnesses

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <span>
#include <string>
#include <string_view>

namespace primerpicker::fuzz {

// The fuzzer's bytes, as a span (the harness entry point gets a
// pointer and a size)
inline auto bytes(std::uint8_t const* const data, std::size_t const size) noexcept
    -> std::span<std::uint8_t const> {
    return {data, size};
}

// The fuzzer's bytes as text (a copy)
inline auto text(std::span<std::uint8_t const> const input) -> std::string {
    auto characters = std::string{};
    characters.reserve(input.size());
    std::ranges::transform(input, std::back_inserter(characters),
                           [](std::uint8_t const byte) { return static_cast<char>(byte); });
    return characters;
}

// Stop on a broken invariant: libFuzzer records the input as a crash
inline auto check(bool const invariant) noexcept -> void {
    if (not invariant) {
        std::abort();
    }
}

// Upper-case letters only (what the sequence readers keep)
inline auto is_upper_case_letters(std::string_view const sequence) noexcept -> bool {
    return std::ranges::all_of(sequence, [](char const letter) { return letter >= 'A' and letter <= 'Z'; });
}

}  // namespace primerpicker::fuzz
