#include "core/analysis/legacy_verdict.hpp"
#include "core/report/legacy_report.hpp"
#include "core/search/legacy_find_pair.hpp"

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cstddef>
#include <string>

namespace {

constexpr auto date = std::chrono::year_month_day{std::chrono::year{2026}, std::chrono::October, std::chrono::day{20}};

auto report_for(primerpicker::search::FindPairSettings const& settings, std::size_t const molecule_size)
    -> std::string {
    return primerpicker::report::legacy_find_pair_report(
        {.molecule_name = "TLIN", .molecule_size = molecule_size, .settings = settings, .best_pairs = {}},
        primerpicker::analysis::LegacyCriteria{}, date);
}

}  // namespace

// reports/find_pair/TLIN_15_560.PRN (circular/TLIN.GB, 600 bases, targets
// 15 and 560, offsets 100): the report prints the offsets cut at the
// molecule's ends (tests/legacy/README.md, "Circular molecules")
TEST_CASE("Find Pair report: offsets cut at the molecule's ends", "[report][legacy]") {
    auto const report = report_for({.target_a = 15, .target_b = 560}, 600);
    CHECK(report.find("   Primer A (normal strand): Target bp:    15;   Start: Target - 14 bases\r\n") !=
          std::string::npos);
    CHECK(report.find("   Primer B (compl. strand): Target bp:   560;   End:   Target + 40 bases\r\n") !=
          std::string::npos);
}

// reports/FINDPAIR.PRN (pairs/TEST600.GB, targets 150 and 400): offsets
// inside the molecule are printed as they are
TEST_CASE("Find Pair report: offsets inside the molecule", "[report][legacy]") {
    auto const report = report_for({.target_a = 150, .target_b = 400}, 600);
    CHECK(report.find("Start: Target - 100 bases\r\n") != std::string::npos);
    CHECK(report.find("End:   Target + 100 bases\r\n") != std::string::npos);
}
