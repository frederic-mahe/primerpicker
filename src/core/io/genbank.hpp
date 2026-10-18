#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace primerpicker::io {

struct Molecule {
    std::string name;
    std::string sequence;  // uppercase letters
};

enum struct GenbankError : std::uint8_t { none, missing_locus, missing_origin };

// C++23 refactoring: replace with std::expected<Molecule, GenbankError>
struct GenbankResult {
    std::optional<Molecule> molecule;
    GenbankError error = GenbankError::none;
};

namespace detail {

// Splits a text into lines, without their line ends (LF or CRLF)
class Lines {
public:
    constexpr explicit Lines(std::string_view const text) noexcept : rest_{text} {}

    [[nodiscard]] constexpr auto next() noexcept -> std::optional<std::string_view> {
        if (rest_.empty()) {
            return std::nullopt;
        }
        auto const end = rest_.find('\n');
        auto line = rest_;
        if (end == std::string_view::npos) {
            rest_ = {};
        } else {
            line.remove_suffix(line.size() - end);
            rest_.remove_prefix(end + 1);
        }
        if (line.ends_with('\r')) {
            line.remove_suffix(1);
        }
        return line;
    }

private:
    std::string_view rest_;
};

// the first word after a keyword such as LOCUS
constexpr auto first_word_after(std::string_view line, std::string_view const keyword) noexcept
    -> std::string_view {
    line.remove_prefix(keyword.size());
    auto const start = line.find_first_not_of(" \t");
    if (start == std::string_view::npos) {
        return {};
    }
    line.remove_prefix(start);
    auto const end = line.find_first_of(" \t");
    if (end != std::string_view::npos) {
        line.remove_suffix(line.size() - end);
    }
    return line;
}

constexpr auto to_upper(char const letter) noexcept -> char {
    constexpr auto case_offset = 'a' - 'A';
    return letter >= 'a' and letter <= 'z' ? static_cast<char>(letter - case_offset) : letter;
}

constexpr auto is_letter(char const character) noexcept -> bool {
    return (character >= 'a' and character <= 'z') or (character >= 'A' and character <= 'Z');
}

// letters of an ORIGIN line ("        1 acgtacgtac gt"); digits and
// spaces are dropped
constexpr auto append_letters(std::string& sequence, std::string_view const line) -> void {
    for (auto const character : line) {
        if (is_letter(character)) {
            sequence.push_back(to_upper(character));
        }
    }
}

}  // namespace detail

// Reads the name (first word after LOCUS) and the sequence (letters of
// the ORIGIN section, up to "//") of a GenBank flat file, as written by
// Primer Designer and Clone Manager. Other sections are ignored.
// Not noexcept: the strings allocate.
[[nodiscard]] constexpr auto parse_genbank(std::string_view const text) -> GenbankResult {
    constexpr auto locus = std::string_view{"LOCUS"};
    constexpr auto origin = std::string_view{"ORIGIN"};
    constexpr auto end_of_entry = std::string_view{"//"};
    auto lines = detail::Lines{text};
    auto name = std::optional<std::string>{};
    auto in_origin = false;
    auto has_origin = false;
    auto sequence = std::string{};
    for (auto line = lines.next(); line.has_value(); line = lines.next()) {
        if (in_origin) {
            if (line->starts_with(end_of_entry)) {
                break;
            }
            detail::append_letters(sequence, *line);
            continue;
        }
        if (not name.has_value() and line->starts_with(locus)) {
            name = std::string{detail::first_word_after(*line, locus)};
            continue;
        }
        if (name.has_value() and line->starts_with(origin)) {
            in_origin = true;
            has_origin = true;
        }
    }
    if (not name.has_value()) {
        return {.molecule = std::nullopt, .error = GenbankError::missing_locus};
    }
    if (not has_origin) {
        return {.molecule = std::nullopt, .error = GenbankError::missing_origin};
    }
    return {.molecule = Molecule{.name = std::move(*name), .sequence = std::move(sequence)},
            .error = GenbankError::none};
}

}  // namespace primerpicker::io
