#include "core/io/molecule_file.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/report/legacy_report.hpp"
#include "core/search/legacy_match.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
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
#include <string_view>
#include <utility>
#include <vector>

namespace {

using primerpicker::io::LegacySystemDefaults;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::Strand;

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

auto test600() -> std::vector<Nucleotide> {
    auto const text = read_file(legacy_dir / "pairs" / "TEST600.GB");
    auto const parsed = primerpicker::io::parse_molecule_file({.text = text, .name = "TEST600.GB"});
    auto sequence = parsed.molecule.has_value() ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                                                : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot parse TEST600.GB"};
    }
    return std::move(*sequence);
}

auto match_all_settings() -> LegacySystemDefaults {
    auto const text = read_file(legacy_dir / "match_variations" / "PD2_MATCH_ALL.CFG");
    auto bytes = std::vector<std::byte>(text.size());
    std::ranges::transform(text, bytes.begin(), [](char const character) { return static_cast<std::byte>(character); });
    auto const parsed = primerpicker::io::parse_pd2_cfg(bytes);
    if (not parsed.defaults.has_value()) {
        throw std::runtime_error{"cannot parse PD2_MATCH_ALL.CFG"};
    }
    return *parsed.defaults;
}

struct Recorded {
    std::string_view report;  // in tests/legacy/complement
    std::size_t position_a;   // complement strand
    std::size_t target;
    Strand strand_b;
    bool match_all;  // PD2_MATCH_ALL.CFG, or the built-in settings
};

// the report of Match Active Primer for a complement-strand primer A
// of 20 bases, printed on 10/22/26
auto report_for(std::vector<Nucleotide> const& molecule, Recorded const& recorded) -> std::string {
    auto const defaults = recorded.match_all ? match_all_settings() : LegacySystemDefaults{};
    auto const primer_a = primerpicker::sequence::primer_at(
        molecule, {.position = recorded.position_a, .length = defaults.primer_length, .strand = Strand::complement});
    if (not primer_a.has_value()) {
        throw std::runtime_error{"no primer A"};
    }
    auto const settings = primerpicker::search::MatchSettings{.target = recorded.target,
                                                              .start_offset = defaults.match_start_offset,
                                                              .end_offset = defaults.match_end_offset,
                                                              .strand = recorded.strand_b,
                                                              .length = defaults.primer_length,
                                                              .gc_tolerance_percent = defaults.gc_tolerance_percent,
                                                              .tm_tolerance_celsius = defaults.tm_tolerance_celsius};
    auto const list = primerpicker::search::legacy_match_active_primer(
        {.molecule = molecule, .primer_a = *primer_a}, settings, defaults.criteria);
    constexpr auto printed = std::chrono::year_month_day{std::chrono::year{2026}, std::chrono::October,
                                                         std::chrono::day{22}};
    return primerpicker::report::legacy_match_report({.molecule_name = "TEST600",
                                                      .molecule_size = molecule.size(),
                                                      .primer_a = *primer_a,
                                                      .primer_a_position = recorded.position_a,
                                                      .primer_a_strand = Strand::complement,
                                                      .settings = settings,
                                                      .criteria = defaults.criteria,
                                                      .list = list},
                                                     printed);
}

}  // namespace

// Reference: PRIMER.exe, Match Active Primer with primer A on the
// complement strand of TEST600, printed to disk on 10/22/26
// (tests/legacy/complement/*.PRN, recorded 2026-10-22): primers B on
// the top strand (the form's default for such a primer A) or on the
// complement strand ("Complement strand (Y/N)" answered Y). Primer A's
// row has no '!' (400(C): GC 45, 460(C): GC 50, the top-strand primer at
// 460 having 45)
TEST_CASE("Match Active Primer with a complement-strand primer A: every recorded report", "[search][report][legacy]") {
    auto const molecule = test600();
    for (auto const& recorded :
         {Recorded{.report = "match_c400.PRN", .position_a = 400, .target = 100, .strand_b = Strand::top, .match_all = false},
          Recorded{.report = "match_c400_all.PRN", .position_a = 400, .target = 100, .strand_b = Strand::top, .match_all = true},
          Recorded{.report = "match_c460_all.PRN", .position_a = 460, .target = 300, .strand_b = Strand::top, .match_all = true},
          Recorded{.report = "match_c460.PRN", .position_a = 460, .target = 300, .strand_b = Strand::top, .match_all = false},
          Recorded{.report = "match_c460_200.PRN", .position_a = 460, .target = 200, .strand_b = Strand::top, .match_all = false},
          Recorded{.report = "match_c400_compB.PRN",
                   .position_a = 400,
                   .target = 450,
                   .strand_b = Strand::complement,
                   .match_all = false}}) {
        CAPTURE(recorded.report);
        CHECK(report_for(molecule, recorded) == read_file(legacy_dir / "complement" / recorded.report));
    }
}
