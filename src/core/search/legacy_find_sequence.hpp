#pragma once

#include "core/enzymes/site_search.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "core/sequence/topology.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

// Find Sequence of Primer Designer 2.0 (Molecule-Ops -> View/Edit ->
// Sequence -> Search; tests/legacy/view_edit/): the sites of a sequence
// of IUPAC letters on both strands of the molecule, and their list.

namespace primerpicker::search {

// The sites of LETTERS (IUPAC, either case) on both strands, by
// position, as the restriction sites (enzymes::find_sites: the leftmost
// base on the top strand; a palindrome once); none when LETTERS is empty
// or holds another letter. Not noexcept: vectors allocate.
[[nodiscard]] constexpr auto legacy_find_sequence(std::span<sequence::Nucleotide const> const molecule,
                                                  std::string_view const letters,
                                                  sequence::Topology const topology)
    -> std::optional<std::vector<enzymes::Site>> {
    auto const pattern = enzymes::parse_site(letters);
    if (not pattern.has_value()) {
        return std::nullopt;
    }
    return enzymes::find_sites(molecule, *pattern, enzymes::Mismatches::none, topology);
}

// The list of sites as PRIMER.exe shows it from column 2 of the screen
// (as the other View/Edit views): each site its position right-aligned
// in 10 columns (9 for a line's first site, whose field starts in column
// 1) and C on the complement strand (a space otherwise), six sites a
// line, separated by commas, without trailing spaces ("        5C,        10 ,       113C,
// ..."). Not noexcept: strings allocate.
[[nodiscard]] constexpr auto legacy_site_list_lines(std::span<enzymes::Site const> const sites)
    -> std::vector<std::string> {
    constexpr auto sites_per_line = std::size_t{6};
    constexpr auto position_width = std::size_t{10};
    constexpr auto decimal_base = std::size_t{10};
    auto lines = std::vector<std::string>{};
    for (auto index = std::size_t{0}; index < sites.size(); ++index) {
        if (index % sites_per_line == 0) {
            lines.emplace_back();
        } else {
            lines.back() += ',';
        }
        auto number = std::string{};
        for (auto value = sites[index].position; value > 0; value /= decimal_base) {
            number.insert(number.begin(), static_cast<char>('0' + (value % decimal_base)));
        }
        auto& line = lines.back();
        auto const width = index % sites_per_line == 0 ? position_width - 1 : position_width;
        if (number.size() < width) {
            line.append(width - number.size(), ' ');
        }
        line += number;
        line += sites[index].strand == sequence::Strand::complement ? 'C' : ' ';
    }
    for (auto& line : lines) {
        line.erase(line.find_last_not_of(' ') + 1);
    }
    return lines;
}

}  // namespace primerpicker::search
