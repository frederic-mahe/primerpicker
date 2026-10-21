#include "core/io/molecule_file.hpp"
#include "core/report/legacy_report.hpp"
#include "core/search/legacy_match.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using primerpicker::sequence::Nucleotide;

auto const legacy_dir = std::filesystem::path{PRIMERPICKER_LEGACY_DIR};

auto read_file(std::filesystem::path const& path) -> std::string {
    auto file = std::ifstream{path, std::ios::binary};
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    auto text = std::ostringstream{};
    text << file.rdbuf();
    return std::move(text).str();
}

auto trun() -> std::vector<Nucleotide> {
    auto const text = read_file(legacy_dir / "search" / "TRUN.GB");
    auto const parsed = primerpicker::io::parse_molecule_file({.text = text, .name = "TRUN.GB"});
    auto sequence = parsed.molecule.has_value() ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                                                : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot parse TRUN.GB"};
    }
    return std::move(*sequence);
}

}  // namespace

// Reference: PRIMER.exe, Match Active Primer for the primer at 50 of
// search/TRUN.GB, target 600, "Complement strand (Y/N)" answered N
// (primers B on the top strand), built-in settings: the report printed
// to disk on 10/20/26 (tests/legacy/match_variations/
// TRUN_50_600_top_strand.PRN, recorded 2026-10-20): "Target bp pos:"
// without "(Complementary strand)", rows 552 to 609
TEST_CASE("Match Active Primer with primers B on the top strand", "[search][report][legacy]") {
    auto const molecule = trun();
    constexpr auto position_a = std::size_t{50};
    constexpr auto length = std::size_t{20};
    auto const primer_a = std::span{molecule}.subspan(position_a - 1, length);
    auto const settings = primerpicker::search::MatchSettings{.target = 600, .strand = primerpicker::sequence::Strand::top};
    auto const list = primerpicker::search::legacy_match_active_primer({.molecule = molecule, .primer_a = primer_a},
                                                                       settings, {});
    REQUIRE_FALSE(list.empty());
    CHECK(list.front().position == 552);
    CHECK(list.back().position == 609);
    constexpr auto printed = std::chrono::year_month_day{std::chrono::year{2026}, std::chrono::October,
                                                         std::chrono::day{20}};
    auto const report = primerpicker::report::legacy_match_report({.molecule_name = "TRUN",
                                                                   .molecule_size = molecule.size(),
                                                                   .primer_a = primer_a,
                                                                   .primer_a_position = position_a,
                                                                   .settings = settings,
                                                                   .criteria = {},
                                                                   .list = list},
                                                                  printed);
    CHECK(report == read_file(legacy_dir / "match_variations" / "TRUN_50_600_top_strand.PRN"));
}
