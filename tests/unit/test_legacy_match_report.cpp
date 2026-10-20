#include "core/io/molecule_file.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/report/legacy_report.hpp"
#include "core/search/legacy_match.hpp"
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
using primerpicker::report::legacy_match_report;
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

struct RecordedReport {
    std::string file;
    std::string molecule;  // relative to tests/legacy
    std::size_t primer_a_position;  // top strand
    std::size_t target;
    std::string config;
};

// reports/match_reports.tsv: columns file, molecule, primer_a_position,
// target, config
auto recorded_reports() -> std::vector<RecordedReport> {
    constexpr auto column_count = 5U;
    auto reports = std::vector<RecordedReport>{};
    auto file = std::ifstream{reports_dir / "match_reports.tsv"};
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
                           .molecule = fields[1],
                           .primer_a_position = static_cast<std::size_t>(number(fields[2])),
                           .target = static_cast<std::size_t>(number(fields[3])),
                           .config = fields[4]});
    }
    return reports;
}

auto report_for(RecordedReport const& recorded, std::chrono::year_month_day const& date) -> std::string {
    auto const path = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / recorded.molecule;
    auto const text = read_file(path);
    auto const parsed = primerpicker::io::parse_molecule_file({.text = text, .name = path.filename().string()});
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"cannot parse " + recorded.molecule};
    }
    auto const molecule = nucleotides(parsed.molecule->sequence);
    auto const defaults = defaults_from(recorded.config);
    auto const primer_a = std::span{molecule}.subspan(recorded.primer_a_position - 1, defaults.primer_length);
    auto const settings = primerpicker::search::MatchSettings{.target = recorded.target,
                                                              .start_offset = defaults.match_start_offset,
                                                              .end_offset = defaults.match_end_offset,
                                                              .length = defaults.primer_length,
                                                              .gc_tolerance_percent = defaults.gc_tolerance_percent,
                                                              .tm_tolerance_celsius = defaults.tm_tolerance_celsius};
    auto const list = primerpicker::search::legacy_match_active_primer({.molecule = molecule, .primer_a = primer_a},
                                                                       settings, defaults.criteria);
    return legacy_match_report({.molecule_name = parsed.molecule->name,
                                .molecule_size = molecule.size(),
                                .primer_a = primer_a,
                                .primer_a_position = recorded.primer_a_position,
                                .settings = settings,
                                .criteria = defaults.criteria,
                                .list = list},
                               date);
}

}  // namespace

// Sweep: every Match Active Primer report recorded with Print to Disk
// (tests/legacy/README.md, "The Match Active Primer report"), with the
// date it was printed
TEST_CASE("Match Active Primer reports: every recorded report, byte for byte", "[report][legacy]") {
    auto const reports = recorded_reports();
    REQUIRE(reports.size() >= 4);
    for (auto const& recorded : reports) {
        INFO(recorded.file);
        auto const expected = read_file(reports_dir / recorded.file);
        CHECK(report_for(recorded, report_date(expected)) == expected);
    }
}
