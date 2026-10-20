#pragma once

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_dimer_drawing.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>

namespace primerpicker::analysis {

// The two duplexes that Primer Designer 2.0 prints in its reports
// ("LARGEST DUPLEXES REPORTED IN PRIMER DIMERS"; tests/legacy/README.md,
// "Duplexes of the reports"). The first primer is drawn 5' to 3', the
// second 3' to 5' below it, each base under its partner. Offsets are in
// columns from the first primer's 5' base, and negative when left of it.
struct DuplexPlacement {
    std::ptrdiff_t second_3prime_offset;  // the second primer's 3' base
    std::ptrdiff_t marks_offset;          // the first pair of the duplex
    std::size_t count;                    // pairs in the duplex (':' marks)
};

// "Primer-dimers with 3'end involvement": the duplex of the dimer page
// (legacy_dimer_drawing). Both primers need at least two nucleotides.
[[nodiscard]] constexpr auto legacy_3prime_duplex(PrimerPair const& pair_of_primers) noexcept -> DuplexPlacement {
    assert(pair_of_primers.first.size() >= 2 and pair_of_primers.second.size() >= 2);
    auto const column = static_cast<std::ptrdiff_t>(detail::dimer_drawing_column);
    auto const [second_3prime_column, marks_column, count] = detail::placement(pair_of_primers);
    return {.second_3prime_offset = second_3prime_column - column,
            .marks_offset = marks_column - column,
            .count = count};
}

namespace detail {

// The longest run of pairs on one diagonal (first position i faces
// second position diagonal - i), and where it starts; the first run
// found when several are equally long
struct DiagonalRun {
    std::size_t start;
    std::size_t count;
};

constexpr auto first_longest_run_on_diagonal(PrimerPair const& pair_of_primers, std::size_t const diagonal) noexcept
    -> DiagonalRun {
    auto const& [first, second] = pair_of_primers;
    auto const first_position = diagonal >= second.size() ? diagonal - second.size() + 1 : 0;
    auto const last_position = std::min(diagonal, first.size() - 1);
    auto longest = DiagonalRun{.start = first_position, .count = 0};
    auto current = std::size_t{0};
    for (auto position = first_position; position <= last_position; ++position) {
        current = pair(first[position], second[diagonal - position]) ? current + 1 : 0;
        if (current > longest.count) {
            longest = DiagonalRun{.start = position + 1 - current, .count = current};
        }
    }
    return longest;
}

}  // namespace detail

// "Other Primer-Dimers": the longest run of consecutive pairs in any
// antiparallel alignment (its length is longest_complementary_run);
// among equally long runs, the first one from the alignment of both 5'
// ends on, then from the first primer's 5' end. Both primers must be
// non-empty.
[[nodiscard]] constexpr auto legacy_longest_duplex(PrimerPair const& pair_of_primers) noexcept -> DuplexPlacement {
    assert(not pair_of_primers.first.empty() and not pair_of_primers.second.empty());
    auto const diagonals = pair_of_primers.first.size() + pair_of_primers.second.size() - 1;
    auto best_diagonal = std::size_t{0};
    auto best = detail::first_longest_run_on_diagonal(pair_of_primers, best_diagonal);
    for (auto diagonal = std::size_t{1}; diagonal < diagonals; ++diagonal) {
        auto const run = detail::first_longest_run_on_diagonal(pair_of_primers, diagonal);
        if (run.count > best.count) {
            best = run;
            best_diagonal = diagonal;
        }
    }
    // second position j faces first position diagonal - j: its 3' base
    // (j = size - 1) is that many columns from the first primer's 5' base
    auto const second_last = static_cast<std::ptrdiff_t>(pair_of_primers.second.size()) - 1;
    return {.second_3prime_offset = static_cast<std::ptrdiff_t>(best_diagonal) - second_last,
            .marks_offset = static_cast<std::ptrdiff_t>(best.start),
            .count = best.count};
}

}  // namespace primerpicker::analysis
