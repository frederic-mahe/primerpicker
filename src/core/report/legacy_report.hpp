#pragma once

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_homology.hpp"
#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_pair.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/search/legacy_match.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <algorithm>
#include <cassert>
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

// A primer edited on the Mutagenize Primer screen: the primer made from
// the molecule (5' to 3', top strand) and the mutated primer, as long
struct MutatedPrimerReport {
    std::span<sequence::Nucleotide const> original;
    std::span<sequence::Nucleotide const> mutated;
    std::string_view molecule_name;
    std::size_t molecule_size;
    std::size_t position;  // the primer's 5' end on the molecule
};

// The One Primer report of a mutated primer (I on the Mutagenize Primer
// screen, then Print to Disk): the source molecule, "Mutated", the
// position, then, on a second page, "PRIMER MUTAGENESIS INFORMATION":
// the bases changed and the translations of both primers in three
// frames (tests/legacy/README.md, "Mutated primers"). Both primers need
// at least two nucleotides. Not noexcept: as above.
// C++26 refactoring: make constexpr (std::format)
[[nodiscard]] auto legacy_mutated_inspect_report(MutatedPrimerReport const& mutated,
                                                 analysis::LegacyCriteria const& criteria,
                                                 std::chrono::year_month_day const& date) -> std::string;

// "Primer Design Report:  Primer Pair", the report of the Primer Pair
// Analysis screen (two pages, a form feed between them). Both primers
// need at least two nucleotides. Not noexcept: as above.
// C++26 refactoring: make constexpr (std::format)
[[nodiscard]] auto legacy_pair_report(analysis::PrimerPair const& primers, analysis::LegacyCriteria const& criteria,
                                      std::chrono::year_month_day const& date) -> std::string;
// The same for primers with a 5' extension: the sequence written with a
// dash after the extension, the values of the whole primer (the Tm
// without the extension), the duplexes of the whole primers
// (tests/legacy/extension/pair_ext_report.PRN). Both primers need at
// least two nucleotides after their extension.
[[nodiscard]] auto legacy_pair_report(analysis::LegacyPrimerPair const& extended,
                                      analysis::LegacyCriteria const& criteria,
                                      std::chrono::year_month_day const& date) -> std::string;

// "Primer Design Report:  Multiple Primers", the report of the Multi
// Primer Analysis screen. Each primer needs at least two nucleotides.
// Not noexcept: as above.
// C++26 refactoring: make constexpr (std::format)
[[nodiscard]] auto legacy_multi_primer_report(analysis::MultiPrimerFields const& fields,
                                              analysis::LegacyCriteria const& criteria,
                                              std::chrono::year_month_day const& date) -> std::string;
// The same for primers with a 5' extension, as the Primer Pair report
// (tests/legacy/extension/multi_ext_report.PRN)
[[nodiscard]] auto legacy_multi_primer_report(analysis::LegacyMultiPrimerFields const& fields,
                                              analysis::LegacyCriteria const& criteria,
                                              std::chrono::year_month_day const& date) -> std::string;

// What the Find Pair report shows: the molecule, the search, and the
// pairs printed, best first. PRIMER.exe asks "N Pairs found.  Print
// best: [M]", M being N up to legacy_report_best_pairs.
constexpr auto legacy_report_best_pairs = std::size_t{20};
// The field takes 4 digits (tests/legacy/README.md, "Reports of more
// than 20 pairs")
constexpr auto legacy_report_max_best_pairs = std::size_t{9999};

// The answer PRIMER.exe proposes for PAIR_COUNT pairs found
[[nodiscard]] constexpr auto legacy_print_best_default(std::size_t const pair_count) noexcept -> std::size_t {
    return std::min(pair_count, legacy_report_best_pairs);
}

// The pairs printed for the answer BEST (0 to 9999): the best BEST of the
// RANKED pairs, all of them when there are fewer; none for 0 (PRIMER.exe
// then prints nothing)
[[nodiscard]] constexpr auto legacy_printed_pairs(std::span<search::PairCandidate const> const ranked,
                                                  std::size_t const best) noexcept
    -> std::span<search::PairCandidate const> {
    assert(best <= legacy_report_max_best_pairs);
    return ranked.first(std::min(ranked.size(), best));
}

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

// What the Uniqueness and False Priming reports show: the primer (taken
// from the molecule, top strand, at position), the molecule, and the
// sites in the order of the screen
struct HomologyReport {
    std::string_view molecule_name;
    std::size_t molecule_size;
    std::span<sequence::Nucleotide const> primer;
    std::size_t position;  // 1-based, the primer's 5' end
    analysis::HomologyKind kind;
    analysis::HomologyTmFilter tm_filter;
    std::span<analysis::HomologySite const> sites;
};

// "Uniqueness of Primer Report" or "False Priming Report", with the
// text map of the molecule (tests/legacy/README.md, "The Uniqueness and
// False Priming reports"). The molecule must not be empty. Not
// noexcept: as above.
// C++26 refactoring: make constexpr (std::format)
[[nodiscard]] auto legacy_homology_report(HomologyReport const& homology, std::chrono::year_month_day const& date)
    -> std::string;

// What the Match Active Primer report shows: the molecule, primer A
// (top strand, at its position), the search, and the list of primers B
struct MatchReport {
    std::string_view molecule_name;
    std::size_t molecule_size;
    std::span<sequence::Nucleotide const> primer_a;
    std::size_t primer_a_position;
    sequence::Strand primer_a_strand = sequence::Strand::top;
    search::MatchSettings settings;
    analysis::LegacyCriteria criteria;  // for primer A's '!' flags
    std::span<search::MatchCandidate const> list;
};

// "Primer Design Report:  Match Active Primer", the report of the Find
// Matching Primer list (tests/legacy/README.md, "The Match Active Primer
// report"). Primer A needs at least two nucleotides. Not noexcept: as
// above.
// C++26 refactoring: make constexpr (std::format)
[[nodiscard]] auto legacy_match_report(MatchReport const& match, std::chrono::year_month_day const& date)
    -> std::string;

}  // namespace primerpicker::report
