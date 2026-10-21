#include "core/analysis/legacy_verdict.hpp"
#include "core/report/legacy_report.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

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
using primerpicker::report::legacy_mutated_inspect_report;
using primerpicker::report::MutatedPrimerReport;
using primerpicker::sequence::Nucleotide;

auto const recordings_dir = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "mutagenize_enter";

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

auto report_for(std::string_view const mutated_letters, std::chrono::year_month_day const& date) -> std::string {
    auto const original = nucleotides("CTCAAAGGCTCAGACATCAG");
    auto const mutated = nucleotides(mutated_letters);
    constexpr auto position = std::size_t{51};
    constexpr auto molecule_size = std::size_t{600};
    return legacy_mutated_inspect_report(MutatedPrimerReport{.original = original,
                                                             .mutated = mutated,
                                                             .molecule_name = "TEST600",
                                                             .molecule_size = molecule_size,
                                                             .position = position},
                                         LegacyCriteria{}, date);
}

}  // namespace

// Reference: PRIMER.exe, Modify -> Mutagenize on the primer at 51 of
// pairs/TEST600.GB (CTCAAAGGCTCAGACATCAG), GAATTC typed over its first
// six bases, then I and Print to Disk, built-in criteria
// (tests/legacy/mutagenize_enter/TEST600_51_GAATTC_i_report.PRN, recorded
// 2026-10-21)
TEST_CASE("the One Primer report of a mutated primer, six bases changed", "[report][legacy][mutagenize]") {
    auto const expected = read_file(recordings_dir / "TEST600_51_GAATTC_i_report.PRN");
    CHECK(report_for("GAATTCGGCTCAGACATCAG", report_date(expected)) == expected);
}

// Same reference, the first base changed to G, then S (saved) and I
// (TEST600_51_G_s_i_report.PRN), and with F1, d, Enter before I
// (TEST600_51_G_f1_d_enter_i_report.PRN: the same report)
TEST_CASE("the One Primer report of a mutated primer, one base changed", "[report][legacy][mutagenize]") {
    for (auto const* const name : {"TEST600_51_G_s_i_report.PRN", "TEST600_51_G_f1_d_enter_i_report.PRN"}) {
        INFO(name);
        auto const expected = read_file(recordings_dir / name);
        CHECK(report_for("GTCAAAGGCTCAGACATCAG", report_date(expected)) == expected);
    }
}
