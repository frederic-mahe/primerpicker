#pragma once

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_homology.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <algorithm>
#include <cstddef>
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
    std::size_t start_offset_a = 100;
    std::size_t target_b = 0;
    std::size_t end_offset_b = 100;
    std::size_t length = 20;  // bases, without 5' extension here
    int gc_tolerance_percent = 5;
    int tm_tolerance_celsius = 10;
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

struct FindPairResult {
    std::size_t possible_a;  // primers A that meet every criterion
    std::size_t possible_b;
    std::vector<PairCandidate> pairs;  // by position of A, then of B
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

// The primers of a window that meet every criterion
inline auto possible_primers(std::span<sequence::Nucleotide const> const molecule, std::size_t const first,
                             std::size_t const last, sequence::Strand const strand, std::size_t const length,
                             analysis::LegacyCriteria const& criteria) -> std::vector<Possible> {
    auto possible = std::vector<Possible>{};
    for (auto position = first; position <= last; ++position) {
        auto primer = sequence::primer_at(molecule, {.position = position, .length = length, .strand = strand});
        if (not primer.has_value()) {
            continue;
        }
        auto const actual = analysis::legacy_actual_values(*primer);
        if (not analysis::evaluate(actual, criteria).all()) {
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
// fewer A:B 3' end matches than the criterion. Not sorted by score
// yet: how false priming enters the score is not determined. Not
// constexpr: the false-priming search uses std::ranges::stable_sort;
// not noexcept: vectors allocate.
[[nodiscard]] inline auto legacy_find_pair(std::span<sequence::Nucleotide const> const molecule,
                                           FindPairSettings const& settings, analysis::LegacyCriteria const& criteria)
    -> FindPairResult {
    auto const first_a = settings.target_a > settings.start_offset_a ? settings.target_a - settings.start_offset_a : 1;
    auto const primers_a = detail::possible_primers(molecule, first_a, settings.target_a, sequence::Strand::top,
                                                    settings.length, criteria);
    auto const primers_b = detail::possible_primers(molecule, settings.target_b, settings.target_b + settings.end_offset_b,
                                                    sequence::Strand::complement, settings.length, criteria);
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

}  // namespace primerpicker::search
