#include "core/analysis/legacy_homology.hpp"
#include "core/io/molecule_file.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/report/legacy_report.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

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
using primerpicker::report::legacy_homology_report;
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
    std::size_t position;  // of the primer's 5' end, top strand
    primerpicker::analysis::HomologyKind kind;
    std::string config;
};

// reports/homology_reports_copies.tsv: columns file, molecule, position, kind
// (uniqueness or false_priming), config, note
auto recorded_reports() -> std::vector<RecordedReport> {
    constexpr auto column_count = 5U;
    auto reports = std::vector<RecordedReport>{};
    auto file = std::ifstream{reports_dir / "homology_reports_copies.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        if (fields.size() != column_count or (fields[3] != "uniqueness" and fields[3] != "false_priming")) {
            throw std::runtime_error{"unexpected line: " + line};
        }
        reports.push_back({.file = fields[0],
                           .molecule = fields[1],
                           .position = static_cast<std::size_t>(number(fields[2])),
                           .kind = fields[3] == "uniqueness" ? primerpicker::analysis::HomologyKind::uniqueness
                                                            : primerpicker::analysis::HomologyKind::false_priming,
                           .config = fields[4]});
    }
    return reports;
}

auto report_for(RecordedReport const& recorded, std::chrono::year_month_day const& date) -> std::string {
    using primerpicker::analysis::HomologyOrder;
    using primerpicker::analysis::HomologyTmFilter;
    auto const path = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / recorded.molecule;
    auto const text = read_file(path);
    auto const name = path.filename().string();
    auto const parsed = primerpicker::io::parse_molecule_file({.text = text, .name = name});
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"cannot parse " + recorded.molecule};
    }
    auto const molecule = nucleotides(parsed.molecule->sequence);
    auto const defaults = defaults_from(recorded.config);
    auto const primer = std::span{molecule}.subspan(recorded.position - 1, defaults.primer_length);
    auto const settings = primerpicker::analysis::HomologySettings{
        .tm_filter = defaults.exclude_tm_at_most_zero ? HomologyTmFilter::exclude_at_most_zero : HomologyTmFilter::keep_all,
        .order = defaults.sort_homologies_by_tm ? HomologyOrder::by_tm : HomologyOrder::by_position};
    auto const own_site = primerpicker::sequence::PrimerLocation{
        .position = recorded.position, .length = primer.size(), .strand = primerpicker::sequence::Strand::top};
    auto const sites = primerpicker::analysis::legacy_homology_sites(primer, molecule, own_site, recorded.kind, settings);
    return legacy_homology_report({.molecule_name = parsed.molecule->name,
                                   .molecule_size = molecule.size(),
                                   .primer = primer,
                                   .position = recorded.position,
                                   .kind = recorded.kind,
                                   .tm_filter = settings.tm_filter,
                                   .sites = sites},
                                  date);
}

}  // namespace

// Sweep (tests/legacy/README.md, "Copies of the primer elsewhere"):
// copies of the primer elsewhere in the molecule, full-length ones
// included, are listed; only the primer's own site is left out
TEST_CASE("Uniqueness and False Priming reports: copies of the primer elsewhere", "[report][legacy]") {
    auto const reports = recorded_reports();
    REQUIRE(reports.size() >= 6);
    for (auto const& recorded : reports) {
        INFO(recorded.file);
        auto const expected = read_file(reports_dir / recorded.file);
        CHECK(report_for(recorded, report_date(expected)) == expected);
    }
}
