#pragma once

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/search/legacy_find_pair.hpp"

#include <chrono>
#include <cstddef>
#include <span>
#include <string>
#include <string_view>

namespace primerpicker::report {

// The text reports of Primer Designer 2.0 (P -> Print to Disk;
// tests/legacy/README.md, "Reports"), byte for byte: lines end with CR
// LF, and the date is the computer's (MM/DD/YY).
constexpr auto legacy_line_end = std::string_view{"\r\n"};

// "Primer Design Report:  One Primer", the report of the Primer
// Analysis screen. The primer needs at least two nucleotides after its
// extension. Not noexcept: the string allocates, and std::format may
// throw.
// C++26 refactoring: make constexpr (std::format)
[[nodiscard]] auto legacy_inspect_report(analysis::LegacyPrimer const& primer,
                                         analysis::LegacyCriteria const& criteria,
                                         std::chrono::year_month_day const& date) -> std::string;

// "Primer Design Report:  Primer Pair", the report of the Primer Pair
// Analysis screen (two pages, a form feed between them). Both primers
// need at least two nucleotides. Not noexcept: as above.
// C++26 refactoring: make constexpr (std::format)
[[nodiscard]] auto legacy_pair_report(analysis::PrimerPair const& primers, analysis::LegacyCriteria const& criteria,
                                      std::chrono::year_month_day const& date) -> std::string;

// "Primer Design Report:  Multiple Primers", the report of the Multi
// Primer Analysis screen. Each primer needs at least two nucleotides.
// Not noexcept: as above.
// C++26 refactoring: make constexpr (std::format)
[[nodiscard]] auto legacy_multi_primer_report(analysis::MultiPrimerFields const& fields,
                                              analysis::LegacyCriteria const& criteria,
                                              std::chrono::year_month_day const& date) -> std::string;

// What the Find Pair report shows: the molecule, the search, and the
// pairs printed, best first. PRIMER.exe asks "N Pairs found.  Print
// best: [M]", M being N up to legacy_report_best_pairs.
constexpr auto legacy_report_best_pairs = std::size_t{20};

struct FindPairReport {
    std::string_view molecule_name;
    std::size_t molecule_size;
    search::FindPairSettings settings;
    std::span<search::PairCandidate const> best_pairs;
};

// "Primer Design Report:  Find Pair", the report of the Find Primer
// Pair list. Not noexcept: as above.
// C++26 refactoring: make constexpr (std::format)
[[nodiscard]] auto legacy_find_pair_report(FindPairReport const& find_pair, analysis::LegacyCriteria const& criteria,
                                           std::chrono::year_month_day const& date) -> std::string;

}  // namespace primerpicker::report
