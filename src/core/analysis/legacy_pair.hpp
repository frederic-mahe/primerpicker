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

}  // namespace primerpicker::analysis
