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
using primerpicker::analysis::LegacyPrimer;
using primerpicker::report::legacy_inspect_report;
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
    std::string primer;
    std::string extension;
    std::string config;
    std::string duplexes;
};

// reports/inspect_reports.tsv: columns file, primer, extension, config,
// duplexes ("-" for none)
auto recorded_reports() -> std::vector<RecordedReport> {
    constexpr auto column_count = 5U;
    auto reports = std::vector<RecordedReport>{};
    auto file = std::ifstream{reports_dir / "inspect_reports.tsv"};
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
        reports.push_back({.file = fields[0],
                           .primer = fields[1],
                           .extension = fields[2] == "-" ? std::string{} : fields[2],
                           .config = fields[3],
                           .duplexes = fields[4]});
    }
    return reports;
}

auto report_for(RecordedReport const& recorded, std::chrono::year_month_day const& date) -> std::string {
    auto sequence = nucleotides(recorded.extension);
    auto const primer = nucleotides(recorded.primer);
    sequence.insert(sequence.end(), primer.begin(), primer.end());
    return legacy_inspect_report(LegacyPrimer{.sequence = sequence, .extension_length = recorded.extension.size()},
                                 criteria_from(recorded.config), date);
}

}  // namespace

// Sweep: every Inspect report recorded with Print to Disk
// (tests/legacy/README.md, "Reports"), with the date it was printed
TEST_CASE("Inspect reports: every recorded report, byte for byte", "[report][legacy]") {
    auto const reports = recorded_reports();
    REQUIRE(reports.size() >= 30);
    for (auto const& recorded : reports) {
        if (recorded.duplexes == "no-pairs") {
            continue;  // a duplex without pairs: see the test below
        }
        INFO(recorded.file);
        auto const expected = read_file(reports_dir / recorded.file);
        CHECK(report_for(recorded, report_date(expected)) == expected);
    }
}

TEST_CASE("Inspect reports: the first line and the duplexes of INSPECT.PRN", "[report][legacy]") {
    // reports/INSPECT.PRN (Direct Entry GTAAAACGACGGCCAGT, printed on
    // 10/19/26)
    auto const sequence = nucleotides("GTAAAACGACGGCCAGT");
    auto const report = legacy_inspect_report(LegacyPrimer{.sequence = sequence, .extension_length = 0},
                                              LegacyCriteria{},
                                              std::chrono::year_month_day{std::chrono::year{2026},
                                                                          std::chrono::October, std::chrono::day{19}});
    CHECK(report.starts_with("Primer Design Report:  One Primer                     Report Date: 10/19/26\r\n"));
    CHECK(report.ends_with("     Other Primer-Dimers:\r\n"
                           "\r\n"
                           "          5'-GTAAAACGACGGCCAGT-3'\r\n"
                           "                       ::::\r\n"
                           "                 3'-TGACCGGCAGCAAAATG-5'\r\n"));
}
