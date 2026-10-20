#pragma once

#include "core/analysis/legacy_builtin_defaults.hpp"
#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_homology.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace primerpicker::search {

// The "Define Primer Pair Constraints" form of Primer Designer 2.0.
// Defaults: the built-in system defaults (tests/legacy/README.md,
// "Built-in defaults"). Primer A is on the top strand, searched from
// target_a - start_offset_a to target_a; primer B on the complement
// strand, from target_b to target_b + end_offset_b; positions are 5'
// ends (see sequence::PrimerLocation).
struct FindPairSettings {
    std::size_t target_a = 0;
    std::size_t start_offset_a = analysis::builtin::find_pair_start_offset;
    std::size_t target_b = 0;
    std::size_t end_offset_b = analysis::builtin::find_pair_end_offset;
    std::size_t length = analysis::builtin::primer_length;  // bases, without 5' extension here
    int gc_tolerance_percent = analysis::builtin::gc_tolerance_percent;
    int tm_tolerance_celsius = analysis::builtin::tm_tolerance_celsius;
};

struct PairCandidate {
    std::size_t position_a;
    std::size_t position_b;
    analysis::LegacyActualValues primer_a;
    analysis::LegacyActualValues primer_b;
    std::size_t max_dimer_3prime_matches;  // largest of A:A, B:B and A:B
    std::size_t max_dimer_adjacent_bases;
    std::optional<int> false_priming_tm_a;  // highest false-priming site Tm, if any
    std::optional<int> false_priming_tm_b;
};

// The limits of Find Pair in Primer Designer 2.0 (tests/legacy/README.md,
// "Limits on primers and pairs (recorded 2026-10-19)"): a window with
// this many possible primers stops the search, and the pair count stops
// at this many pairs
constexpr auto legacy_max_possible_primers = std::size_t{100};
constexpr auto legacy_max_pairs = std::size_t{1000};

enum struct FindPairLimit : std::uint8_t {
    none,
    // a window has legacy_max_possible_primers possible primers: no
    // pair is searched ("More than 100 sites found.  Suggest you
    // redefine matching limits.")
    too_many_primers,
    // more pairs match than legacy_max_pairs: the first ones are kept
    // ("More than 1000 possible pairs found!"; PRIMER.exe hangs instead)
    pairs_capped,
};

struct FindPairResult {
    std::size_t possible_a;  // primers A that meet every criterion
    std::size_t possible_b;  // 0 when primer A reached the limit
    std::vector<PairCandidate> pairs;  // by position of A, then of B
    FindPairLimit limit = FindPairLimit::none;
};

namespace detail {

struct Possible {
    std::size_t position;
    std::vector<sequence::Nucleotide> sequence;
    analysis::LegacyActualValues actual;
    std::optional<int> false_priming_tm;
};

constexpr auto difference(int const lhs, int const rhs) noexcept -> int {
    return lhs > rhs ? lhs - rhs : rhs - lhs;
}

// The run of one base at the 3' end of a primer (1 when the last two
// bases differ). Requires a non-empty primer.
constexpr auto run_at_3prime_end(std::span<sequence::Nucleotide const> const primer) noexcept -> std::size_t {
    assert(not primer.empty());
    auto const last = primer.back();
    auto const run = std::ranges::find_if(primer.rbegin(), primer.rend(),
                                          [last](sequence::Nucleotide const nucleotide) { return nucleotide != last; });
    return static_cast<std::size_t>(std::distance(primer.rbegin(), run));
}

// The primers of a window that meet every criterion; the count stops at
// legacy_max_possible_primers
inline auto possible_primers(std::span<sequence::Nucleotide const> const molecule, std::size_t const first,
                             std::size_t const last, sequence::Strand const strand, std::size_t const length,
                             analysis::LegacyCriteria const& criteria) -> std::vector<Possible> {
    auto possible = std::vector<Possible>{};
    for (auto position = first; position <= last and possible.size() < legacy_max_possible_primers; ++position) {
        auto primer = sequence::primer_at(molecule, {.position = position, .length = length, .strand = strand});
        if (not primer.has_value()) {
            continue;
        }
        auto const actual = analysis::legacy_actual_values(*primer);
        // unlike the Primer Analysis screen, Find Pair also rejects a run
        // at the 3' end that reaches the runs criterion
        // (tests/legacy/README.md, "Find Pair and runs at the 3' end")
        if (not analysis::evaluate(actual, criteria).all() or run_at_3prime_end(*primer) >= criteria.reject_run_length) {
            continue;
        }
        auto const sites = analysis::legacy_homology_sites(*primer, molecule, analysis::HomologyKind::false_priming);
        auto const highest_tm = sites.empty() ? std::nullopt : std::optional{sites.front().tm_celsius};
        possible.push_back({.position = position, .sequence = std::move(*primer), .actual = actual, .false_priming_tm = highest_tm});
    }
    return possible;
}

}  // namespace detail

// Find Pair of Primer Designer 2.0 (tests/legacy/README.md, "Find
// Pair"): the possible primers A and B (every criterion met), and the
// pairs whose GC % and Tm match within the tolerances and that have
// fewer A:B 3' end matches than the criterion, in generation order;
// legacy_rank_pairs sorts them as the Find Pair list. The search stops
// at the limits of FindPairLimit. Not
// constexpr: the false-priming search uses std::ranges::stable_sort;
// not noexcept: vectors allocate.
[[nodiscard]] inline auto legacy_find_pair(std::span<sequence::Nucleotide const> const molecule,
                                           FindPairSettings const& settings, analysis::LegacyCriteria const& criteria)
    -> FindPairResult {
    auto const first_a = settings.target_a > settings.start_offset_a ? settings.target_a - settings.start_offset_a : 1;
    auto const last_a = sequence::window_end(settings.target_a, 0, molecule.size());
    auto const primers_a =
        detail::possible_primers(molecule, first_a, last_a, sequence::Strand::top, settings.length, criteria);
    if (primers_a.size() == legacy_max_possible_primers) {
        return {.possible_a = primers_a.size(), .possible_b = 0, .pairs = {}, .limit = FindPairLimit::too_many_primers};
    }
    auto const last_b = sequence::window_end(settings.target_b, settings.end_offset_b, molecule.size());
    auto const primers_b = detail::possible_primers(molecule, settings.target_b, last_b,
                                                    sequence::Strand::complement, settings.length, criteria);
    if (primers_b.size() == legacy_max_possible_primers) {
        return {.possible_a = primers_a.size(),
                .possible_b = primers_b.size(),
                .pairs = {},
                .limit = FindPairLimit::too_many_primers};
    }
    auto result = FindPairResult{.possible_a = primers_a.size(), .possible_b = primers_b.size(), .pairs = {}};
    for (auto const& primer_a : primers_a) {
        for (auto const& primer_b : primers_b) {
            auto const pair = analysis::PrimerPair{.first = primer_a.sequence, .second = primer_b.sequence};
            auto const inter_3prime = analysis::dimer_at_3prime_end(pair);
            auto const matches =
                detail::difference(primer_a.actual.gc_percent, primer_b.actual.gc_percent) <= settings.gc_tolerance_percent and
                detail::difference(primer_a.actual.tm_celsius, primer_b.actual.tm_celsius) <= settings.tm_tolerance_celsius and
                inter_3prime < criteria.reject_dimer_3prime_matches;
            if (not matches) {
                continue;
            }
            if (result.pairs.size() == legacy_max_pairs) {
                result.limit = FindPairLimit::pairs_capped;
                return result;
            }
            auto const& actual_a = primer_a.actual;
            auto const& actual_b = primer_b.actual;
            result.pairs.push_back(PairCandidate{
                .position_a = primer_a.position,
                .position_b = primer_b.position,
                .primer_a = actual_a,
                .primer_b = actual_b,
                .max_dimer_3prime_matches = std::max({actual_a.dimer_3prime_matches, actual_b.dimer_3prime_matches, inter_3prime}),
                .max_dimer_adjacent_bases = std::max({actual_a.dimer_adjacent_bases, actual_b.dimer_adjacent_bases,
                                                      analysis::longest_complementary_run(pair)}),
                .false_priming_tm_a = primer_a.false_priming_tm,
                .false_priming_tm_b = primer_b.false_priming_tm});
        }
    }
    return result;
}

// The "Find-Pair Weighting Factors" (F4 in the pair constraints form):
// 0 ignores a term, 1 to 5 weigh it. Defaults: the built-in system
// defaults (tests/legacy/README.md, "Built-in defaults").
struct FindPairWeights {
    int gc_mismatch = analysis::builtin::find_pair_weight;
    int tm_mismatch = analysis::builtin::find_pair_weight;
    int dimer_3prime = analysis::builtin::find_pair_weight;
    int dimer_any = analysis::builtin::find_pair_weight;
    int false_priming = analysis::builtin::find_pair_weight;
};

namespace detail {

constexpr auto max_weight = 5;

constexpr auto is_weight(int const weight) noexcept -> bool {
    return weight >= 0 and weight <= max_weight;
}

}  // namespace detail

// Score of a pair in Find Pair of Primer Designer 2.0, lower is better
// (tests/legacy/README.md, "Find Pair"): the sum of weight x term, the
// terms being the GC % difference, the Tm difference, the max 3' dimer
// matches, the max adjacent bases, and the sum of the false-priming
// site Tm of primers A and B (0 for a primer without site)
[[nodiscard]] constexpr auto legacy_pair_score(PairCandidate const& pair, FindPairWeights const& weights) noexcept
    -> int {
    assert(detail::is_weight(weights.gc_mismatch) and detail::is_weight(weights.tm_mismatch) and
           detail::is_weight(weights.dimer_3prime) and detail::is_weight(weights.dimer_any) and
           detail::is_weight(weights.false_priming));
    auto const gc_mismatch = detail::difference(pair.primer_a.gc_percent, pair.primer_b.gc_percent);
    auto const tm_mismatch = detail::difference(pair.primer_a.tm_celsius, pair.primer_b.tm_celsius);
    auto const dimer_3prime = static_cast<int>(pair.max_dimer_3prime_matches);
    auto const dimer_any = static_cast<int>(pair.max_dimer_adjacent_bases);
    auto const false_priming = pair.false_priming_tm_a.value_or(0) + pair.false_priming_tm_b.value_or(0);
    return (weights.gc_mismatch * gc_mismatch) + (weights.tm_mismatch * tm_mismatch) +
           (weights.dimer_3prime * dimer_3prime) + (weights.dimer_any * dimer_any) +
           (weights.false_priming * false_priming);
}

// The Find Pair list of Primer Designer 2.0: the pairs of
// legacy_find_pair (in generation order: by position of A, then of B),
// by increasing score; pairs of equal score keep their generation order.
// C++26 refactoring: make constexpr (std::ranges::stable_sort)
[[nodiscard]] inline auto legacy_rank_pairs(std::vector<PairCandidate> pairs, FindPairWeights const& weights)
    -> std::vector<PairCandidate> {
    auto const score = [&weights](PairCandidate const& pair) { return legacy_pair_score(pair, weights); };
    std::ranges::stable_sort(pairs, std::less{}, score);
    return pairs;
}

}  // namespace primerpicker::search
