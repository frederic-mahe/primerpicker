#pragma once

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_verdict.hpp"

#include <cstddef>
#include <cstdlib>

namespace primerpicker::analysis {

// The pair-specific values of the Primer Pair Analysis screen of
// Primer Designer 2.0 (tests/legacy/README.md, "Primer pair analysis"):
// balance differences (shown as "N base difference", "N % difference",
// "N°C difference" when not zero) and A:B dimers
struct LegacyPairValues {
    int length_difference;
    int gc_difference_percent;
    int tm_difference_celsius;
    std::size_t dimer_3prime_matches;
    std::size_t dimer_adjacent_bases;
};

// Both primers need at least two nucleotides (without 5' extension)
[[nodiscard]] constexpr auto legacy_pair_values(PrimerPair const& primers) noexcept -> LegacyPairValues {
    auto const first = legacy_actual_values(primers.first);
    auto const second = legacy_actual_values(primers.second);
    auto const length_difference = static_cast<std::ptrdiff_t>(primers.first.size()) -
                                   static_cast<std::ptrdiff_t>(primers.second.size());
    auto const absolute = [](auto const value) { return value < 0 ? -value : value; };
    return {.length_difference = static_cast<int>(absolute(length_difference)),
            .gc_difference_percent = absolute(first.gc_percent - second.gc_percent),
            .tm_difference_celsius = absolute(first.tm_celsius - second.tm_celsius),
            .dimer_3prime_matches = dimer_at_3prime_end(primers),
            .dimer_adjacent_bases = longest_complementary_run(primers)};
}

// The "!" flags of one primer on the Primer Pair Analysis screen
struct LegacyPairCautions {
    bool gc;             // GC % outside the criteria range
    bool tm;             // Tm outside the criteria range
    bool gc_run_3prime;  // "Runs of G/C at 3' end"
    bool runs;           // "Runs of bases (any)"
    bool hairpin;        // "Hairpins"
};

// Flags as recorded on 42 pair screens: GC % and Tm as on the Primer
// Analysis screen, a hairpin as there; but runs and 3' G/C runs are
// flagged only above the "reject if >=" value (a run of 3 is not
// flagged with the default 3, while the Primer Analysis screen fails
// it)
[[nodiscard]] constexpr auto legacy_pair_cautions(LegacyActualValues const& actual,
                                                  LegacyCriteria const& criteria) noexcept -> LegacyPairCautions {
    auto const verdict = evaluate(actual, criteria);
    return {.gc = not verdict.gc,
            .tm = not verdict.tm,
            .gc_run_3prime = actual.gc_run_3prime > criteria.reject_gc_run_3prime,
            .runs = actual.longest_run > criteria.reject_run_length,
            .hairpin = not verdict.no_hairpin};
}

struct LegacyDimerFlags {
    bool three_prime;  // "3' end matches"
    bool adjacent;     // "Any adj homol bases"
};

// Dimer flags (A:A, B:B or A:B): the criteria's "reject if >=" values
[[nodiscard]] constexpr auto legacy_dimer_flags(std::size_t const three_prime_matches,
                                                std::size_t const adjacent_bases,
                                                LegacyCriteria const& criteria) noexcept -> LegacyDimerFlags {
    return {.three_prime = three_prime_matches >= criteria.reject_dimer_3prime_matches,
            .adjacent = adjacent_bases >= criteria.reject_dimer_adjacent_bases};
}

}  // namespace primerpicker::analysis
