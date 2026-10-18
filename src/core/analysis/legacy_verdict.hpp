#pragma once

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_hairpin.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/thermo/legacy_tm.hpp"
#include "core/thermo/units.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <optional>
#include <span>

namespace primerpicker::analysis {

// Single-primer criteria of Primer Designer 2.0 ("Criteria" screen).
// Default values: the built-in defaults of PRIMER.exe
// (tests/legacy/README.md, "Built-in defaults").
struct LegacyCriteria {
    int gc_min_percent = 50;
    int gc_max_percent = 60;
    int tm_min_celsius = 55;
    int tm_max_celsius = 80;
    // a hairpin with an energy below the cutoff fails the criterion
    thermo::KilocaloriesPerMole hairpin_cutoff = thermo::KilocaloriesPerMole{0.0};
    std::size_t reject_dimer_3prime_matches = 3;  // reject if >=
    std::size_t reject_dimer_adjacent_bases = 7;  // reject if >=
    std::size_t reject_run_length = 3;            // reject if >=
    std::size_t reject_gc_run_3prime = 3;         // reject if >=
};

// The "Actual" column of the Primer Analysis screen
struct LegacyActualValues {
    int gc_percent;
    int tm_celsius;  // 0 when the Tm is not positive ("?" elsewhere)
    std::optional<thermo::KilocaloriesPerMole> hairpin;  // the most stable hairpin, if any
    std::size_t dimer_3prime_matches;
    std::size_t dimer_adjacent_bases;
    std::size_t longest_run;
    std::size_t gc_run_3prime;
};

// The "Meets Criteria" column; all() is "MEETS CRITERIA SELECTED"
struct LegacyVerdict {
    bool gc;
    bool tm;
    bool no_hairpin;
    bool no_3prime_dimer;
    bool no_dimer;
    bool no_runs;
    bool no_3prime_gc_run;

    [[nodiscard]] constexpr auto all() const noexcept -> bool {
        return gc and tm and no_hairpin and no_3prime_dimer and no_dimer and no_runs and no_3prime_gc_run;
    }
};

// GC content as shown: a truncated percentage. Requires a non-empty
// primer.
[[nodiscard]] constexpr auto legacy_gc_percent(std::span<sequence::Nucleotide const> const primer) noexcept
    -> int {
    assert(not primer.empty());
    constexpr auto percent = 100;
    auto const gc_count = std::ranges::count_if(primer, [](sequence::Nucleotide const nucleotide) {
        return nucleotide == sequence::Nucleotide::cytosine or nucleotide == sequence::Nucleotide::guanine;
    });
    return static_cast<int>((gc_count * percent) / std::ssize(primer));
}

// The "Actual" column for a primer (without 5' extension), computed
// as Primer Designer does. Requires at least two nucleotides.
[[nodiscard]] constexpr auto legacy_actual_values(std::span<sequence::Nucleotide const> const primer) noexcept
    -> LegacyActualValues {
    assert(primer.size() >= 2);
    auto const hairpin = legacy_hairpin(primer);
    // the Tm is truncated, and shown as 0 when not positive
    auto const primer_tm = thermo::legacy_primer_tm(primer).degrees;
    auto const truncated_tm = static_cast<int>(primer_tm);  // toward zero: the floor when positive
    return {.gc_percent = legacy_gc_percent(primer),
            .tm_celsius = primer_tm > 0.0 ? truncated_tm : 0,
            .hairpin = hairpin.has_value() ? std::optional{hairpin->energy} : std::nullopt,
            .dimer_3prime_matches = self_dimer_at_3prime_end(primer),
            .dimer_adjacent_bases = longest_self_complementary_run(primer),
            .longest_run = longest_run_before_3prime_end(primer),
            .gc_run_3prime = gc_run_at_3prime_end(primer)};
}

// Range limits are inclusive (PRIMER.exe shows YES at 50 and 60 %, 55
// and 80 °C, the default limits)
[[nodiscard]] constexpr auto evaluate(LegacyActualValues const& actual, LegacyCriteria const& criteria) noexcept
    -> LegacyVerdict {
    auto const within = [](int const value, int const minimum, int const maximum) {
        return value >= minimum and value <= maximum;
    };
    auto const hairpin_fails = actual.hairpin.has_value() and
                               actual.hairpin->kcal_per_mol < criteria.hairpin_cutoff.kcal_per_mol;
    return {.gc = within(actual.gc_percent, criteria.gc_min_percent, criteria.gc_max_percent),
            .tm = within(actual.tm_celsius, criteria.tm_min_celsius, criteria.tm_max_celsius),
            .no_hairpin = not hairpin_fails,
            .no_3prime_dimer = actual.dimer_3prime_matches < criteria.reject_dimer_3prime_matches,
            .no_dimer = actual.dimer_adjacent_bases < criteria.reject_dimer_adjacent_bases,
            .no_runs = actual.longest_run < criteria.reject_run_length,
            .no_3prime_gc_run = actual.gc_run_3prime < criteria.reject_gc_run_3prime};
}

}  // namespace primerpicker::analysis
