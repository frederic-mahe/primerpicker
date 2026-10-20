#pragma once

#include "core/analysis/legacy_builtin_defaults.hpp"
#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <cstddef>
#include <span>
#include <utility>
#include <vector>

namespace primerpicker::search {

// The "Define Matching Primer" form of Primer Designer 2.0. Defaults:
// the built-in system defaults (tests/legacy/README.md, "Built-in
// defaults") and the form's own (complement strand).
struct MatchSettings {
    std::size_t target = 0;           // 1-based molecule position
    std::size_t start_offset = analysis::builtin::match_start_offset;  // search from target - start_offset
    std::size_t end_offset = analysis::builtin::match_end_offset;      // to target + end_offset
    sequence::Strand strand = sequence::Strand::complement;
    std::size_t length = analysis::builtin::primer_length;           // bases, without 5' extension here
    int gc_tolerance_percent = analysis::builtin::gc_tolerance_percent;  // match GC % of primer A within
    int tm_tolerance_celsius = analysis::builtin::tm_tolerance_celsius;  // match Tm of primer A within
};

// The molecule to search, and the active primer (primer A)
struct MatchInput {
    std::span<sequence::Nucleotide const> molecule;
    std::span<sequence::Nucleotide const> primer_a;
};

struct MatchCandidate {
    std::size_t position;  // 5' end, as in PrimerLocation
    std::size_t distance;  // "bps from target"
    std::vector<sequence::Nucleotide> sequence;
    analysis::LegacyActualValues actual;
    std::size_t inter_dimer_3prime_matches;  // with primer A
    std::size_t inter_dimer_adjacent_bases;
};

namespace detail {

constexpr auto absolute_difference(int const lhs, int const rhs) noexcept -> int {
    return lhs > rhs ? lhs - rhs : rhs - lhs;
}

}  // namespace detail

// Match Active Primer of Primer Designer 2.0 (tests/legacy/README.md,
// "Match Active Primer"): the primers B of the search window that meet
// every criterion, whose GC % and Tm are within the tolerances of
// primer A's, and that have fewer A:B 3' end matches than the
// criterion; sorted by position. Not noexcept: vectors allocate.
[[nodiscard]] constexpr auto legacy_match_active_primer(MatchInput const& input, MatchSettings const& settings,
                                                        analysis::LegacyCriteria const& criteria)
    -> std::vector<MatchCandidate> {
    auto const& [molecule, primer_a] = input;
    auto const actual_a = analysis::legacy_actual_values(primer_a);
    auto const first = settings.target > settings.start_offset ? settings.target - settings.start_offset : 1;
    auto const last = sequence::window_end(settings.target, settings.end_offset, molecule.size());
    auto candidates = std::vector<MatchCandidate>{};
    for (auto position = first; position <= last; ++position) {
        auto primer_b = sequence::primer_at(
            molecule, {.position = position, .length = settings.length, .strand = settings.strand});
        if (not primer_b.has_value()) {
            continue;
        }
        auto const actual_b = analysis::legacy_actual_values(*primer_b);
        auto const pair = analysis::PrimerPair{.first = primer_a, .second = *primer_b};
        auto const inter_3prime = analysis::dimer_at_3prime_end(pair);
        // as Find Pair, also a run at the 3' end that reaches the runs
        // criterion (tests/legacy/README.md, "Find Pair and runs at the 3'
        // end")
        auto const matches = analysis::evaluate(actual_b, criteria).all() and
                             analysis::run_at_3prime_end(*primer_b) < criteria.reject_run_length and
                             detail::absolute_difference(actual_a.gc_percent, actual_b.gc_percent) <=
                                 settings.gc_tolerance_percent and
                             detail::absolute_difference(actual_a.tm_celsius, actual_b.tm_celsius) <=
                                 settings.tm_tolerance_celsius and
                             inter_3prime < criteria.reject_dimer_3prime_matches;
        if (not matches) {
            continue;
        }
        auto const distance = position > settings.target ? position - settings.target : settings.target - position;
        candidates.push_back(MatchCandidate{.position = position,
                                            .distance = distance,
                                            .sequence = std::move(*primer_b),
                                            .actual = actual_b,
                                            .inter_dimer_3prime_matches = inter_3prime,
                                            .inter_dimer_adjacent_bases = analysis::longest_complementary_run(pair)});
    }
    return candidates;
}

}  // namespace primerpicker::search
