#pragma once

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_verdict.hpp"

#include <chrono>
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

}  // namespace primerpicker::report
