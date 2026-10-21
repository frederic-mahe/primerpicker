#pragma once

#include "core/sequence/nucleotide.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <span>
#include <string>
#include <string_view>
#include <vector>

// The views of Molecule-Ops -> View/Edit of Primer Designer 2.0
// (tests/legacy/view_edit/): the text of each view from column 2 of the
// screen, one string per line.

namespace primerpicker::io {

namespace detail {

constexpr auto view_bases_per_line = std::size_t{60};
constexpr auto view_bases_per_block = std::size_t{10};
constexpr auto view_number_width = std::size_t{7};

// the lines of a file, without their line ends (LF or CR LF)
constexpr auto file_lines(std::string_view text) -> std::vector<std::string_view> {
    auto lines = std::vector<std::string_view>{};
    while (not text.empty()) {
        auto const end = text.find('\n');
        auto line = text.substr(0, end);
        if (line.ends_with('\r')) {
            line.remove_suffix(1);
        }
        lines.push_back(line);
        if (end == std::string_view::npos) {
            break;
        }
        text.remove_prefix(end + 1);
    }
    return lines;
}

}  // namespace detail

// The Sequence view: lines of 60 bases in blocks of 10, after the line's
// first position right-aligned in 7 columns and two spaces
// ("      1  GGTCATTTTA AATGGACAAC ..."). Not noexcept: strings allocate.
[[nodiscard]] constexpr auto legacy_sequence_lines(std::span<sequence::Nucleotide const> const molecule)
    -> std::vector<std::string> {
    auto lines = std::vector<std::string>{};
    for (auto start = std::size_t{0}; start < molecule.size(); start += detail::view_bases_per_line) {
        constexpr auto decimal_base = std::size_t{10};
        auto number = std::string{};
        for (auto value = start + 1; value > 0; value /= decimal_base) {
            number.insert(number.begin(), static_cast<char>('0' + (value % decimal_base)));
        }
        auto line = std::string(detail::view_number_width - std::min(number.size(), detail::view_number_width), ' ');
        line += number;
        line += "  ";
        auto const end = std::min(molecule.size(), start + detail::view_bases_per_line);
        for (auto index = start; index < end; ++index) {
            if (index > start and (index - start) % detail::view_bases_per_block == 0) {
                line += ' ';
            }
            line += sequence::to_letter(molecule[index]);
        }
        lines.push_back(std::move(line));
    }
    return lines;
}

// The Annotations view: the lines of a GenBank file before its
// sequence, ORIGIN included (none for a file without ORIGIN). Not
// noexcept: strings allocate.
[[nodiscard]] constexpr auto legacy_annotation_lines(std::string_view const file_text) -> std::vector<std::string> {
    auto lines = std::vector<std::string>{};
    for (auto const line : detail::file_lines(file_text)) {
        lines.emplace_back(line);
        if (line.starts_with("ORIGIN")) {
            return lines;
        }
    }
    return {};
}

// The Features view: the lines of the FEATURES table of a GenBank file
// (after its FEATURES line, before ORIGIN), as written; "-- no features
// --" when there are none. Not noexcept: strings allocate.
[[nodiscard]] constexpr auto legacy_feature_lines(std::string_view const file_text) -> std::vector<std::string> {
    auto lines = std::vector<std::string>{};
    auto inside = false;
    for (auto const line : detail::file_lines(file_text)) {
        if (line.starts_with("ORIGIN")) {
            break;
        }
        if (inside) {
            lines.emplace_back(line);
        }
        inside = inside or line.starts_with("FEATURES");
    }
    auto const blank = [](std::string const& line) { return line.find_first_not_of(' ') == std::string::npos; };
    if (std::ranges::all_of(lines, blank)) {
        return {"-- no features --"};
    }
    return lines;
}

}  // namespace primerpicker::io
