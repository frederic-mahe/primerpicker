#pragma once

#include "core/sequence/nucleotide.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <ranges>
#include <span>

namespace primerpicker::analysis {

// Single-primer criteria of Primer Designer 2.0, reverse-engineered
// from recorded values (tests/legacy/README.md, "Primer criteria").
// Counts are in nucleotides.

// Two primers facing each other in antiparallel alignments; for a
// self-dimer, both are the same primer
struct PrimerPair {
    std::span<sequence::Nucleotide const> first;
    std::span<sequence::Nucleotide const> second;
};

namespace detail {

using sequence::Nucleotide;

constexpr auto is_gc(Nucleotide const nucleotide) noexcept -> bool {
    return nucleotide == Nucleotide::cytosine or nucleotide == Nucleotide::guanine;
}

constexpr auto pair(Nucleotide const lhs, Nucleotide const rhs) noexcept -> bool {
    return sequence::complement(lhs) == rhs;
}

// First position i faces second position j = diagonal - i. Longest run
// of complementary pairs along one diagonal.
constexpr auto longest_run_on_diagonal(PrimerPair const& pair_of_primers, std::size_t const diagonal) noexcept
    -> std::size_t {
    auto const& [first, second] = pair_of_primers;
    auto const first_position = diagonal >= second.size() ? diagonal - second.size() + 1 : 0;
    auto const last_position = std::min(diagonal, first.size() - 1);
    auto longest = std::size_t{0};
    auto current = std::size_t{0};
    for (auto position = first_position; position <= last_position; ++position) {
        current = pair(first[position], second[diagonal - position]) ? current + 1 : 0;
        longest = std::max(longest, current);
    }
    return longest;
}

// Consecutive complementary pairs from the 3' end of the first primer,
// facing the second primer from position first_partner_position on.
constexpr auto pairs_from_3prime_end(PrimerPair const& pair_of_primers, std::size_t const first_partner_position) noexcept
    -> std::size_t {
    auto const& [primer, partner] = pair_of_primers;
    auto count = std::size_t{0};
    auto position = primer.size() - 1;
    for (auto partner_position = first_partner_position; partner_position < partner.size(); ++partner_position) {
        if (not pair(primer[position], partner[partner_position])) {
            break;
        }
        ++count;
        if (position == 0) {
            break;
        }
        --position;
    }
    return count;
}

// The 3'-end count of the first primer on the second, which must extend
// at least one base beyond that 3' end (as a template would)
constexpr auto directional_dimer_at_3prime_end(PrimerPair const& pair_of_primers) noexcept -> std::size_t {
    auto longest = std::size_t{0};
    for (auto partner_position = std::size_t{1}; partner_position < pair_of_primers.second.size(); ++partner_position) {
        longest = std::max(longest, pairs_from_3prime_end(pair_of_primers, partner_position));
    }
    return longest;
}

}  // namespace detail

// "Dimers, adjacent homologous bases" between two primers: the longest
// run of consecutive complementary pairs, in any antiparallel alignment.
// Both primers must be non-empty.
[[nodiscard]] constexpr auto longest_complementary_run(PrimerPair const& pair_of_primers) noexcept -> std::size_t {
    assert(not pair_of_primers.first.empty() and not pair_of_primers.second.empty());
    auto longest = std::size_t{0};
    auto const diagonals = pair_of_primers.first.size() + pair_of_primers.second.size() - 1;
    for (auto diagonal = std::size_t{0}; diagonal < diagonals; ++diagonal) {
        longest = std::max(longest, detail::longest_run_on_diagonal(pair_of_primers, diagonal));
    }
    return longest;
}

// "3' dimers" between two primers: the larger of the 3'-end counts of
// each primer on the other (see self_dimer_at_3prime_end). Both primers
// must be non-empty.
[[nodiscard]] constexpr auto dimer_at_3prime_end(PrimerPair const& pair_of_primers) noexcept -> std::size_t {
    assert(not pair_of_primers.first.empty() and not pair_of_primers.second.empty());
    return std::max(detail::directional_dimer_at_3prime_end(pair_of_primers),
                    detail::directional_dimer_at_3prime_end(
                        PrimerPair{.first = pair_of_primers.second, .second = pair_of_primers.first}));
}

// "Runs": the longest run of one nucleotide, not counting the run at
// the 3' end (zero for a homopolymer)
[[nodiscard]] constexpr auto longest_run_before_3prime_end(std::span<sequence::Nucleotide const> const primer) noexcept -> std::size_t {
    if (primer.empty()) {
        return 0;
    }
    // start of the run at the 3' end
    auto const last_run = std::ranges::find_if_not(primer | std::views::reverse,
                                                   [&primer](sequence::Nucleotide const nucleotide) {
                                                       return nucleotide == primer.back();
                                                   });
    auto const before = primer.first(static_cast<std::size_t>(std::ranges::distance(last_run, primer.rend())));
    auto longest = std::size_t{0};
    for (auto run = before.begin(); run != before.end();) {
        auto const next = std::ranges::find_if_not(run, before.end(), [run](sequence::Nucleotide const nucleotide) {
            return nucleotide == *run;
        });
        longest = std::max(longest, static_cast<std::size_t>(std::ranges::distance(run, next)));
        run = next;
    }
    return longest;
}

// "3' GC runs": the number of G or C at the 3' end
[[nodiscard]] constexpr auto gc_run_at_3prime_end(std::span<sequence::Nucleotide const> const primer) noexcept -> std::size_t {
    auto const reversed = primer | std::views::reverse;
    return static_cast<std::size_t>(std::ranges::distance(reversed.begin(), std::ranges::find_if_not(reversed, detail::is_gc)));
}

// "Dimers, adjacent homologous bases": the longest run of consecutive
// complementary pairs between the primer and a copy of itself, in any
// antiparallel alignment
[[nodiscard]] constexpr auto longest_self_complementary_run(std::span<sequence::Nucleotide const> const primer) noexcept -> std::size_t {
    if (primer.empty()) {
        return 0;
    }
    return longest_complementary_run(PrimerPair{.first = primer, .second = primer});
}

// "3' dimers, matches at 3' end": the number of consecutive
// complementary pairs from the primer's 3' end, in the best
// antiparallel alignment with a copy of itself where the copy extends
// beyond that 3' end (as a template would)
[[nodiscard]] constexpr auto self_dimer_at_3prime_end(std::span<sequence::Nucleotide const> const primer) noexcept -> std::size_t {
    if (primer.empty()) {
        return 0;
    }
    return detail::directional_dimer_at_3prime_end(PrimerPair{.first = primer, .second = primer});
}

}  // namespace primerpicker::analysis
