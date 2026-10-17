#pragma once

#include "core/sequence/nucleotide.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <ranges>
#include <span>

namespace primerpicker::analysis {

// Single-primer criteria of Primer Designer 2.0, reverse-engineered
// from recorded values (tests/legacy/README.md, "Primer criteria").
// Counts are in nucleotides.

namespace detail {

using sequence::Nucleotide;

constexpr auto is_gc(Nucleotide const nucleotide) noexcept -> bool {
    return nucleotide == Nucleotide::cytosine or nucleotide == Nucleotide::guanine;
}

constexpr auto pair(std::span<Nucleotide const> const primer, std::size_t const lhs, std::size_t const rhs) noexcept
    -> bool {
    return sequence::complement(primer[lhs]) == primer[rhs];
}

// In an antiparallel alignment of the primer with a copy of itself,
// primer position i faces copy position j = diagonal - i.
// Longest run of complementary pairs along one diagonal.
constexpr auto longest_run_on_diagonal(std::span<Nucleotide const> const primer, std::size_t const diagonal) noexcept
    -> std::size_t {
    auto const last = primer.size() - 1;
    auto const first_position = diagonal > last ? diagonal - last : 0;
    auto const last_position = std::min(diagonal, last);
    auto longest = std::size_t{0};
    auto current = std::size_t{0};
    for (auto position = first_position; position <= last_position; ++position) {
        current = pair(primer, position, diagonal - position) ? current + 1 : 0;
        longest = std::max(longest, current);
    }
    return longest;
}

// Consecutive complementary pairs from the primer's 3' end, facing
// copy position first_copy_position and onwards.
constexpr auto pairs_from_3prime_end(std::span<Nucleotide const> const primer, std::size_t const first_copy_position) noexcept
    -> std::size_t {
    auto count = std::size_t{0};
    auto position = primer.size() - 1;
    for (auto copy_position = first_copy_position; copy_position < primer.size(); ++copy_position) {
        if (not pair(primer, position, copy_position)) {
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

}  // namespace detail

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
    auto longest = std::size_t{0};
    auto const diagonals = (2 * primer.size()) - 1;
    for (auto diagonal = std::size_t{0}; diagonal < diagonals; ++diagonal) {
        longest = std::max(longest, detail::longest_run_on_diagonal(primer, diagonal));
    }
    return longest;
}

// "3' dimers, matches at 3' end": the number of consecutive
// complementary pairs from the primer's 3' end, in the best
// antiparallel alignment with a copy of itself where the copy extends
// beyond that 3' end (as a template would)
[[nodiscard]] constexpr auto self_dimer_at_3prime_end(std::span<sequence::Nucleotide const> const primer) noexcept -> std::size_t {
    auto longest = std::size_t{0};
    // the copy must extend beyond the primer's 3' end: the copy position
    // facing that end is 1 or more
    for (auto copy_position = std::size_t{1}; copy_position < primer.size(); ++copy_position) {
        longest = std::max(longest, detail::pairs_from_3prime_end(primer, copy_position));
    }
    return longest;
}

}  // namespace primerpicker::analysis
