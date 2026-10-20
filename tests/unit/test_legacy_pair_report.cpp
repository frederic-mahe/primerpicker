#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/report/legacy_report.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using primerpicker::analysis::LegacyCriteria;
using primerpicker::analysis::PrimerPair;
using primerpicker::report::legacy_pair_report;
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

auto criteria_from(std::string_view const config) -> LegacyCriteria {
    if (config == "-") {
        return LegacyCriteria{};
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
    return result.defaults->criteria;
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

struct RecordedReport {
    std::string file;
    std::string primer_a;
    std::string primer_b;
    std::string config;
    std::string duplexes;
};

// reports/pair_reports.tsv: columns file, primer_a, primer_b, config,
// duplexes
auto recorded_reports() -> std::vector<RecordedReport> {
    constexpr auto column_count = 5U;
    auto reports = std::vector<RecordedReport>{};
    auto file = std::ifstream{reports_dir / "pair_reports.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        if (fields.size() != column_count) {
            throw std::runtime_error{"unexpected line: " + line};
        }
        reports.push_back(
            {.file = fields[0], .primer_a = fields[1], .primer_b = fields[2], .config = fields[3], .duplexes = fields[4]});
    }
    return reports;
}

}  // namespace

// Sweep: every Primer Pair report recorded with Print to Disk
// (tests/legacy/README.md, "Reports"), with the date it was printed
TEST_CASE("Pair reports: every recorded report, byte for byte", "[report][legacy]") {
    auto const reports = recorded_reports();
    REQUIRE(not reports.empty());
    for (auto const& recorded : reports) {
        if (recorded.duplexes == "no-pairs") {
            continue;  // a duplex without pairs: not reproduced (README)
        }
        INFO(recorded.file);
        auto const expected = read_file(reports_dir / recorded.file);
        auto const primer_a = nucleotides(recorded.primer_a);
        auto const primer_b = nucleotides(recorded.primer_b);
        CHECK(legacy_pair_report(PrimerPair{.first = primer_a, .second = primer_b}, criteria_from(recorded.config),
                                 report_date(expected)) == expected);
    }
}

TEST_CASE("Pair reports: CR and form feed before page 2 (PAIR.PRN)", "[report][legacy]") {
    // reports/PAIR.PRN: GTAAAACGACGGCCAGT (A), AACAGCTATGACCATG (B),
    // printed on 10/20/26
    auto const primer_a = nucleotides("GTAAAACGACGGCCAGT");
    auto const primer_b = nucleotides("AACAGCTATGACCATG");
    auto const report = legacy_pair_report(PrimerPair{.first = primer_a, .second = primer_b}, LegacyCriteria{},
                                           std::chrono::year_month_day{std::chrono::year{2026},
                                                                       std::chrono::October, std::chrono::day{20}});
    CHECK(report.find("for duplexes)\r\n\r\n\r\fPrimer Design Report:  Primer Pair  -- Page 2         "
                          "Report Date: 10/20/26\r\n") != std::string::npos);  // C++23 refactoring: use contains()
}
