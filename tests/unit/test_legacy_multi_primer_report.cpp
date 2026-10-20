#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/report/legacy_report.hpp"
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
#include <vector>

namespace {

using primerpicker::analysis::LegacyCriteria;
using primerpicker::analysis::MultiPrimerFields;
using primerpicker::report::legacy_multi_primer_report;
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
    return std::string{std::istreambuf_iterator<char>{file}, std::istreambuf_iterator<char>{}};
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

constexpr auto field_count = primerpicker::analysis::legacy_multi_primer_fields;

struct RecordedReport {
    std::string file;
    std::array<std::string, field_count> fields;  // empty: an empty field
    std::string config;
};

// reports/multi_reports.tsv: columns file, selection (as typed), the
// primers of fields A to D ("-" for an empty field), config
auto recorded_reports() -> std::vector<RecordedReport> {
    constexpr auto column_count = 7U;
    constexpr auto first_field_column = 2U;
    auto reports = std::vector<RecordedReport>{};
    auto file = std::ifstream{reports_dir / "multi_reports.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto columns = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto column = std::string{}; std::getline(stream, column, '\t');) {
            columns.push_back(column == "-" ? std::string{} : column);
        }
        if (columns.size() != column_count) {
            throw std::runtime_error{"unexpected line: " + line};
        }
        auto recorded = RecordedReport{.file = columns[0], .fields = {}, .config = columns.back()};
        for (auto index = std::size_t{0}; index < field_count; ++index) {
            recorded.fields.at(index) = columns.at(first_field_column + index);
        }
        reports.push_back(recorded);
    }
    return reports;
}

}  // namespace

// Sweep: every Multi Primer report recorded with Print to Disk
// (tests/legacy/README.md, "Reports"), with the date it was printed
TEST_CASE("Multi Primer reports: every recorded report, byte for byte", "[report][legacy]") {
    auto const reports = recorded_reports();
    REQUIRE(not reports.empty());
    for (auto const& recorded : reports) {
        INFO(recorded.file);
        auto primers = std::array<std::vector<Nucleotide>, field_count>{};
        auto fields = MultiPrimerFields{};
        for (auto index = std::size_t{0}; index < field_count; ++index) {
            if (recorded.fields.at(index).empty()) {
                continue;
            }
            primers.at(index) = nucleotides(recorded.fields.at(index));
            fields.at(index) = std::span<Nucleotide const>{primers.at(index)};
        }
        auto const config = recorded.config.empty() ? std::string{"-"} : recorded.config;
        auto const expected = read_file(reports_dir / recorded.file);
        CHECK(legacy_multi_primer_report(fields, criteria_from(config), report_date(expected)) == expected);
    }
}
