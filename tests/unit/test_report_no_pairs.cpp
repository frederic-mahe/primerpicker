#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_homology.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/report/legacy_report.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "core/thermo/units.hpp"

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using primerpicker::analysis::LegacyCriteria;
using primerpicker::analysis::LegacyPrimer;
using primerpicker::analysis::PrimerPair;
using primerpicker::report::legacy_inspect_report;
using primerpicker::report::legacy_pair_report;
using primerpicker::sequence::Nucleotide;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

constexpr auto date = std::chrono::year_month_day{std::chrono::year{2026}, std::chrono::October, std::chrono::day{20}};

}  // namespace

// A duplex without pairs. PRIMER.exe writes lines of tens of thousands
// of spaces there (tests/legacy/README.md, "Duplexes without pairs");
// decided with the developer on 2026-10-20: primerpicker draws the
// partner right under the primer, its "3'-" under the "5'-", and an
// empty marks line.
TEST_CASE("Reports: a duplex without pairs is drawn under the primer", "[report]") {
    // AAAAAAAAAAAA pairs with nothing
    auto const primer = nucleotides("AAAAAAAAAAAA");
    auto const report = legacy_inspect_report(LegacyPrimer{.sequence = primer, .extension_length = 0},
                                              LegacyCriteria{}, date);
    CHECK(report.ends_with("     Primer-dimers with 3'end involvement:\r\n"
                           "\r\n"
                           "          5'-AAAAAAAAAAAA-3'\r\n"
                           "\r\n"
                           "          3'-AAAAAAAAAAAA-5'\r\n"
                           "\r\n"
                           "\r\n"
                           "     Other Primer-Dimers:\r\n"
                           "\r\n"
                           "          5'-AAAAAAAAAAAA-3'\r\n"
                           "\r\n"
                           "          3'-AAAAAAAAAAAA-5'\r\n"));
}

TEST_CASE("Reports: only the duplex without pairs changes", "[report]") {
    // ACGCCGTTTTCGGTGT: no pair at the 3' end, a run of 3 elsewhere
    // (reports/inspect/ACGCCGTTTTCGGTGT.PRN, the "Other Primer-Dimers"
    // drawing)
    auto const primer = nucleotides("ACGCCGTTTTCGGTGT");
    auto const report = legacy_inspect_report(LegacyPrimer{.sequence = primer, .extension_length = 0},
                                              LegacyCriteria{}, date);
    CHECK(report.ends_with("     Primer-dimers with 3'end involvement:\r\n"
                           "\r\n"
                           "          5'-ACGCCGTTTTCGGTGT-3'\r\n"
                           "\r\n"
                           "          3'-TGTGGCTTTTGCCGCA-5'\r\n"
                           "\r\n"
                           "\r\n"
                           "     Other Primer-Dimers:\r\n"
                           "\r\n"
                           "          5'-ACGCCGTTTTCGGTGT-3'\r\n"
                           "             :::\r\n"
                           " 3'-TGTGGCTTTTGCCGCA-5'\r\n"));
}

TEST_CASE("Reports: a pair of primers without pairs between them", "[report]") {
    // AAAAAAAAAAAA and CCCCCCCCCCCC (reports/pair/: PRIMER.exe writes
    // spaces there)
    auto const primer_a = nucleotides("AAAAAAAAAAAA");
    auto const primer_b = nucleotides("CCCCCCCCCCCC");
    auto const report = legacy_pair_report(PrimerPair{.first = primer_a, .second = primer_b}, LegacyCriteria{}, date);
    CHECK(report.ends_with("A:B\r\n"
                           "          5'-AAAAAAAAAAAA-3'\r\n"
                           "\r\n"
                           "          3'-CCCCCCCCCCCC-5'\r\n"
                           "\r\n"
                           "\r\n"
                           "---------------------------------------------------------------------------\r\n"));
}

// Found by fuzz_homology_report (2026-10-20): a uniqueness site longer
// than the 19-character bases field pushed the row past the "Bp Pos"
// column, and the padding width wrapped around (a 2 GB allocation). Not
// recorded from PRIMER.exe: primerpicker keeps at least one space before
// the position.
// Since then, PRIMER.exe's layout of sites of 20 and 21 bases was
// recorded (" 5'" in column 1, the bases ending in column 33;
// tests/legacy/README.md, "Copies of the primer elsewhere"): this
// longer site follows it, with at least one space.
TEST_CASE("Reports: a homology site longer than the bases field", "[report]") {
    using primerpicker::analysis::HomologyKind;
    using primerpicker::analysis::HomologySite;
    using primerpicker::analysis::HomologyTmFilter;
    // a 40-base primer, and a site of its first 36 bases at 330
    auto const primer = nucleotides("ACGTACGGTCAGTCCATGCAAGTCGATTGCAGTCCAGTTA");
    auto const sites = std::vector<HomologySite>{
        {.primer_start = 0, .length = 36, .position = 330, .strand = primerpicker::sequence::Strand::top, .tm_celsius = 80}};
    auto const report = primerpicker::report::legacy_homology_report({.molecule_name = "LONG",
                                                                      .molecule_size = 600,
                                                                      .primer = primer,
                                                                      .position = 51,
                                                                      .kind = HomologyKind::uniqueness,
                                                                      .tm_filter = HomologyTmFilter::exclude_at_most_zero,
                                                                      .sites = sites},
                                                                     date);
    CHECK(report.ends_with(" 5' ACGTACGGTCAGTCCATGCAAGTCGATTGCAGTCCA....3'    330            80   \r\n"));
}

// The same wrap-around, in the criteria rows: a settings file may hold
// any hairpin cutoff (a 32-bit float), and a long one pushed the row
// past the "Actual" column. Not recorded from PRIMER.exe: primerpicker
// keeps at least one space before the value.
TEST_CASE("Reports: a criteria setting longer than its column", "[report]") {
    constexpr auto huge_cutoff = 1e38;
    auto criteria = LegacyCriteria{};
    criteria.hairpin_cutoff = primerpicker::thermo::KilocaloriesPerMole{huge_cutoff};
    auto const primer = nucleotides("GTAAAACGACGGCCAGT");
    auto const report =
        legacy_inspect_report(LegacyPrimer{.sequence = primer, .extension_length = 0}, criteria, date);
    auto const row = report.find("No Hairpins      Energy cutoff ");
    REQUIRE(row != std::string::npos);
    auto const row_end = report.find("\r\n", row);
    CHECK(report.substr(row_end - std::string_view{" kcal --           YES"}.size(),
                        std::string_view{" kcal --           YES"}.size()) == " kcal --           YES");
}
