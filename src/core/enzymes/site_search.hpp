#pragma once

#include "core/sequence/iupac.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "core/sequence/topology.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <optional>
#include <ranges>
#include <span>
#include <string_view>
#include <vector>

namespace primerpicker::enzymes {

// A recognition sequence, 5' to 3', one set of nucleotides per letter
using SitePattern = std::vector<sequence::NucleotideSet>;

// The pattern of a recognition sequence written with IUPAC letters (in
// either case); none if it is empty or has other letters. Not
// noexcept: the vector allocates.
[[nodiscard]] constexpr auto parse_site(std::string_view const letters) -> std::optional<SitePattern> {
    if (letters.empty()) {
        return std::nullopt;
    }
    auto pattern = SitePattern{};
    pattern.reserve(letters.size());
    for (auto const letter : letters) {
        auto const set = sequence::from_iupac(letter);
        if (not set.has_value()) {
            return std::nullopt;
        }
        pattern.push_back(*set);
    }
    return pattern;
}

// Primer Designer's "True Sites" (none) and "Almost Sites" (at most one)
enum struct Mismatches : std::uint8_t { none, at_most_one };

// A site on a molecule (tests/legacy/README.md, "True and almost
// sites"): the 1-based top-strand position of its leftmost base, its
// strand (complement: the reverse complement of the recognition
// sequence matches the top strand there), and the offset from the
// leftmost base of the mismatching base, if any
struct Site {
    std::size_t position;
    sequence::Strand strand;
    std::optional<std::size_t> mismatch_offset;

    friend constexpr auto operator==(Site const&, Site const&) -> bool = default;
};

namespace detail {

struct Comparison {
    std::size_t mismatch_count;  // 0, 1 or 2 (2: two or more)
    std::size_t first_mismatch_offset;
};

constexpr auto compare(std::span<sequence::Nucleotide const> const window,
                       std::span<sequence::NucleotideSet const> const pattern) noexcept -> Comparison {
    assert(window.size() == pattern.size());
    constexpr auto enough = std::size_t{2};
    auto comparison = Comparison{.mismatch_count = 0, .first_mismatch_offset = 0};
    for (auto offset = std::size_t{0}; offset < pattern.size() and comparison.mismatch_count < enough; ++offset) {
        if (pattern[offset].contains(window[offset])) {
            continue;
        }
        if (comparison.mismatch_count == 0) {
            comparison.first_mismatch_offset = offset;
        }
        ++comparison.mismatch_count;
    }
    return comparison;
}

// The pattern that the top strand shows where the complement strand
// holds the site
constexpr auto reverse_complement(std::span<sequence::NucleotideSet const> const pattern) -> SitePattern {
    auto reversed = SitePattern{};
    reversed.reserve(pattern.size());
    std::ranges::transform(pattern | std::views::reverse, std::back_inserter(reversed),
                           [](sequence::NucleotideSet const set) { return set.complement(); });
    return reversed;
}

constexpr auto accepted(Comparison const& comparison, Mismatches const mismatches) noexcept -> bool {
    return comparison.mismatch_count == 0 or
           (comparison.mismatch_count == 1 and mismatches == Mismatches::at_most_one);
}

constexpr auto site(std::size_t const offset, sequence::Strand const strand, Comparison const& comparison) noexcept
    -> Site {
    auto const mismatch =
        comparison.mismatch_count == 0 ? std::nullopt : std::optional{comparison.first_mismatch_offset};
    return Site{.position = offset + 1, .strand = strand, .mismatch_offset = mismatch};
}

}  // namespace detail

// The sites of a recognition sequence on both strands of a molecule,
// by position (tests/legacy/README.md, "True and almost sites"); on a
// circular molecule, sites across the origin too. At the same
// position, an exact site comes before a site with a mismatch, then
// top before complement. A palindromic site is found once, on the top
// strand. Not noexcept: vectors allocate.
[[nodiscard]] constexpr auto find_sites(std::span<sequence::Nucleotide const> const molecule,
                                        SitePattern const& pattern, Mismatches const mismatches,
                                        sequence::Topology const topology) -> std::vector<Site> {
    assert(not pattern.empty());
    auto sites = std::vector<Site>{};
    if (molecule.size() < pattern.size()) {
        return sites;
    }
    // a circular molecule is searched on a copy followed by its first
    // bases, from every start position of the molecule
    auto const circular = topology == sequence::Topology::circular;
    auto extended = std::vector<sequence::Nucleotide>(molecule.begin(), molecule.end());
    if (circular) {
        extended.insert(extended.end(), molecule.begin(), std::next(molecule.begin(), std::ssize(pattern) - 1));
    }
    auto const starts = circular ? molecule.size() : molecule.size() - pattern.size() + 1;
    auto const complement_pattern = detail::reverse_complement(pattern);
    auto const palindromic = complement_pattern == pattern;
    for (auto offset = std::size_t{0}; offset < starts; ++offset) {
        auto const window = std::span{extended}.subspan(offset, pattern.size());
        auto const top = detail::compare(window, pattern);
        // a palindrome's complement sites are its top sites: never accepted
        auto const complement = palindromic ? detail::Comparison{.mismatch_count = 2, .first_mismatch_offset = 0}
                                            : detail::compare(window, complement_pattern);
        auto const complement_first = complement.mismatch_count < top.mismatch_count;
        if (complement_first and detail::accepted(complement, mismatches)) {
            sites.push_back(detail::site(offset, sequence::Strand::complement, complement));
        }
        if (detail::accepted(top, mismatches)) {
            sites.push_back(detail::site(offset, sequence::Strand::top, top));
        }
        if (not complement_first and detail::accepted(complement, mismatches)) {
            sites.push_back(detail::site(offset, sequence::Strand::complement, complement));
        }
    }
    return sites;
}

}  // namespace primerpicker::enzymes
