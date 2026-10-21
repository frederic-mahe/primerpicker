#include "core/io/molecule_file.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/report/legacy_report.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using primerpicker::io::LegacySystemDefaults;
using primerpicker::report::legacy_find_pair_report;
using primerpicker::report::legacy_print_best_default;
using primerpicker::report::legacy_printed_pairs;
using primerpicker::sequence::Nucleotide;

auto const reports_dir = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "reports";

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

auto read_file(std::filesystem::path const& path) -> std::string {
    auto file = std::ifstream{path, std::ios::binary};
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    auto text = std::ostringstream{};
    text << file.rdbuf();
    return std::move(text).str();
}

auto defaults_from(std::string_view const config) -> LegacySystemDefaults {
    if (config == "-") {
        return LegacySystemDefaults{};
    }
    auto const text = read_file(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "config" / config);
    auto bytes = std::vector<std::byte>{};
    std::ranges::transform(text, std::back_inserter(bytes), [](char const character) {
        return static_cast<std::byte>(character);
    });
    auto const result = primerpicker::io::parse_pd2_cfg(bytes);
    if (not result.defaults.has_value()) {
        throw std::runtime_error{"cannot parse " + std::string{config}};
    }
    return *result.defaults;
}

auto number(std::string const& text) -> int {
    auto value = 0;
    auto const* const last = std::next(text.data(), std::ssize(text));
    auto const [end, error] = std::from_chars(text.data(), last, value);
    if (error != std::errc{} or end != last) {
        throw std::runtime_error{"not a number: " + text};
    }
    return value;
}

// The date printed on the first line ("Report Date: MM/DD/YY")
auto report_date(std::string_view const report) -> std::chrono::year_month_day {
    constexpr auto label = std::string_view{"Report Date: "};
    constexpr auto field = std::size_t{2};
    constexpr auto twenty_first_century = 2000;
    auto const start = report.find(label);
    if (start == std::string_view::npos) {
        throw std::runtime_error{"no report date"};
    }
    auto const date = std::string{report.substr(start + label.size(), std::string_view{"MM/DD/YY"}.size())};
    return std::chrono::year_month_day{
        std::chrono::year{twenty_first_century + number(date.substr(2 * (field + 1), field))},
        std::chrono::month{static_cast<unsigned>(number(date.substr(0, field)))},
        std::chrono::day{static_cast<unsigned>(number(date.substr(field + 1, field)))}};
}


// A report printed with a typed "Print best: [N]"
struct BestReport {
    std::string_view file;
    std::string_view molecule;
    std::string_view config;
    std::size_t best;
};

// The ranked pairs of Find Pair on MOLECULE, targets 150 and 400, and
// what the report needs
struct Search {
    std::string molecule_name;
    std::size_t molecule_size;
    primerpicker::search::FindPairSettings settings;
    LegacySystemDefaults defaults;
    std::vector<primerpicker::search::PairCandidate> ranked;
};

// a molecule (from tests/legacy) and its settings file ("-": built in)
struct SearchInput {
    std::string_view molecule;
    std::string_view config;
};

auto search(SearchInput const& input) -> Search {
    auto const& [relative_path, config] = input;
    auto const path = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / relative_path;
    auto const name = path.filename().string();
    auto const parsed = primerpicker::io::parse_molecule_file({.text = read_file(path), .name = name});
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"cannot parse " + std::string{relative_path}};
    }
    auto const molecule = nucleotides(parsed.molecule->sequence);
    auto const defaults = defaults_from(config);
    constexpr auto target_a = std::size_t{150};
    constexpr auto target_b = std::size_t{400};
    auto const settings = primerpicker::search::FindPairSettings{.target_a = target_a,
                                                                 .start_offset_a = defaults.find_pair_start_offset,
                                                                 .target_b = target_b,
                                                                 .end_offset_b = defaults.find_pair_end_offset,
                                                                 .length = defaults.primer_length,
                                                                 .gc_tolerance_percent = defaults.gc_tolerance_percent,
                                                                 .tm_tolerance_celsius = defaults.tm_tolerance_celsius};
    auto const result = primerpicker::search::legacy_find_pair(molecule, settings, defaults.criteria);
    return {.molecule_name = parsed.molecule->name,
            .molecule_size = molecule.size(),
            .settings = settings,
            .defaults = defaults,
            .ranked = primerpicker::search::legacy_rank_pairs(result.pairs, defaults.weights)};
}

auto report_for(BestReport const& recorded, std::chrono::year_month_day const& date) -> std::string {
    auto const found = search({.molecule = recorded.molecule, .config = recorded.config});
    auto const& [molecule_name, molecule_size, settings, defaults, ranked] = found;
    return legacy_find_pair_report({.molecule_name = molecule_name,
                                    .molecule_size = molecule_size,
                                    .settings = settings,
                                    .best_pairs = legacy_printed_pairs(ranked, recorded.best)},
                                   defaults.criteria, date);
}

}  // namespace

// Reference: PRIMER.exe, Find Pair on search/TFP.GB (30 pairs) and on
// pairs/TEST600.GB with config/PD2_CRITERIA.CFG (1000 pairs), targets 150
// and 400, Print to Disk with a typed "Print best: [N]"
// (tests/legacy/reports/find_pair_best/, recorded 2026-10-21): the best
// min(N, pairs) rows, one header, no page break
TEST_CASE("Find Pair reports with a typed Print best, byte for byte", "[report][legacy]") {
    constexpr auto recorded = std::array{
        BestReport{.file = "TFP_150_400_best25.PRN", .molecule = "search/TFP.GB", .config = "-", .best = 25},
        BestReport{.file = "TFP_150_400_best30.PRN", .molecule = "search/TFP.GB", .config = "-", .best = 30},
        BestReport{.file = "TFP_150_400_best50.PRN", .molecule = "search/TFP.GB", .config = "-", .best = 50},
        BestReport{.file = "TFP_150_400_best9999.PRN", .molecule = "search/TFP.GB", .config = "-", .best = 9999},
        BestReport{.file = "TEST600_150_400_CRITERIA_best100.PRN",
                   .molecule = "pairs/TEST600.GB",
                   .config = "PD2_CRITERIA.CFG",
                   .best = 100},
        BestReport{.file = "TEST600_150_400_CRITERIA_best1000.PRN",
                   .molecule = "pairs/TEST600.GB",
                   .config = "PD2_CRITERIA.CFG",
                   .best = 1000},
    };
    for (auto const& report : recorded) {
        INFO(report.file);
        auto const expected = read_file(reports_dir / "find_pair_best" / report.file);
        CHECK(report_for(report, report_date(expected)) == expected);
    }
}

// Reference: PRIMER.exe proposes "Print best: [20  ]" for 30 and 1000
// pairs, and "Print best: [M]", M being N up to 20, for fewer
// (tests/legacy/README.md, "Find Pair report" and "Reports of more than
// 20 pairs"); 0 prints nothing (TFP_150_400_best0)
TEST_CASE("Print best: the proposed answer, and 0", "[report][legacy]") {
    CHECK(legacy_print_best_default(30) == 20);
    CHECK(legacy_print_best_default(1000) == 20);
    CHECK(legacy_print_best_default(7) == 7);
    auto const found = search({.molecule = "search/TFP.GB", .config = "-"});
    REQUIRE(found.ranked.size() == 30);
    CHECK(legacy_printed_pairs(found.ranked, 0).empty());
    CHECK(legacy_printed_pairs(found.ranked, 2).size() == 2);
    CHECK(legacy_printed_pairs(found.ranked, 9999).size() == 30);
}
