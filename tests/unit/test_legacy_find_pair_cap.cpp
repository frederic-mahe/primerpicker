#include "core/io/molecule_file.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using primerpicker::io::LegacySystemDefaults;
using primerpicker::search::FindPairLimit;
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

auto settings_file(std::filesystem::path const& path) -> LegacySystemDefaults {
    auto const text = read_file(path);
    auto bytes = std::vector<std::byte>{};
    std::ranges::transform(text, std::back_inserter(bytes),
                           [](char const character) { return static_cast<std::byte>(character); });
    auto const result = primerpicker::io::parse_pd2_cfg(bytes);
    if (not result.defaults.has_value()) {
        throw std::runtime_error{"not a pd2.cfg: " + path.string()};
    }
    return *result.defaults;
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

using Positions = std::pair<std::size_t, std::size_t>;

// the positions A and B of the rows of a Find Pair report, in order
auto report_rows(std::filesystem::path const& path) -> std::vector<Positions> {
    auto lines = std::istringstream{read_file(path)};
    auto rows = std::vector<Positions>{};
    auto in_list = false;
    for (auto line = std::string{}; std::getline(lines, line);) {
        if (line.find("-A-") != std::string::npos) {
            in_list = true;
            continue;
        }
        auto fields = std::istringstream{line};
        auto position_a = std::size_t{0};
        auto position_b = std::size_t{0};
        if (in_list and fields >> position_a >> position_b) {
            rows.emplace_back(position_a, position_b);
        }
    }
    return rows;
}

struct Search {
    std::filesystem::path config;
    std::size_t target_a;
    std::size_t target_b;
};

// the ranked list of primerpicker, positions only
auto ranked_positions(Search const& search) -> std::vector<Positions> {
    auto const defaults = settings_file(search.config);
    auto const settings = primerpicker::search::FindPairSettings{.target_a = search.target_a,
                                                                 .start_offset_a = defaults.find_pair_start_offset,
                                                                 .target_b = search.target_b,
                                                                 .end_offset_b = defaults.find_pair_end_offset,
                                                                 .length = defaults.primer_length,
                                                                 .gc_tolerance_percent = defaults.gc_tolerance_percent,
                                                                 .tm_tolerance_celsius = defaults.tm_tolerance_celsius};
    auto const result = primerpicker::search::legacy_find_pair(test600(), settings, defaults.criteria, defaults.weights);
    CHECK(result.limit == FindPairLimit::pairs_capped);
    CHECK(result.pairs.size() == primerpicker::search::legacy_max_pairs);
    auto positions = std::vector<Positions>{};
    std::ranges::transform(primerpicker::search::legacy_rank_pairs(result.pairs, defaults.weights),
                           std::back_inserter(positions), [](primerpicker::search::PairCandidate const& pair) {
                               return Positions{pair.position_a, pair.position_b};
                           });
    return positions;
}

}  // namespace

// Reference: PRIMER.exe, Find Pair on TEST600.GB, targets 150 and 400,
// config/PD2_CRITERIA.CFG, Print best 1000: the 1000 rows are the best
// 1000 of the first 1001 pairs generated (148/421, the 1001st, is
// listed; 71/422 is not) (tests/legacy/reports/find_pair_best/
// TEST600_150_400_CRITERIA_best1000.PRN, recorded 2026-10-21)
TEST_CASE("Find Pair: the best 1000 of the first 1001 pairs", "[search][legacy]") {
    auto const ranked = ranked_positions({.config = legacy_dir / "config" / "PD2_CRITERIA.CFG", .target_a = 150, .target_b = 400});
    auto const recorded =
        report_rows(legacy_dir / "reports" / "find_pair_best" / "TEST600_150_400_CRITERIA_best1000.PRN");
    REQUIRE(recorded.size() == 1000);
    CHECK(ranked == recorded);
    // C++23 refactoring: use std::ranges::contains
    CHECK(std::ranges::find(ranked, Positions{148, 421}) != ranked.end());
    CHECK(std::ranges::find(ranked, Positions{71, 422}) == ranked.end());
}

// Reference: the same search with the weights 1 0 1 3 1 and targets 260
// and 400, a case designed so that the two models disagree: the 1001st
// pair generated, 260/450, is row 7 of the best 20
// (tests/legacy/find_pair_trun/TEST600_260_400_CRITERIA_W10131.PRN,
// recorded 2026-10-21; its settings file is an input)
TEST_CASE("Find Pair: the 1001st pair can be among the best", "[search][legacy]") {
    auto const ranked = ranked_positions(
        {.config = legacy_dir / "find_pair_trun" / "PD2_CRITERIA_W10131.CFG", .target_a = 260, .target_b = 400});
    auto const recorded = report_rows(legacy_dir / "find_pair_trun" / "TEST600_260_400_CRITERIA_W10131.PRN");
    constexpr auto best = std::size_t{20};
    REQUIRE(recorded.size() == best);
    REQUIRE(ranked.size() >= best);
    CHECK(std::vector<Positions>(ranked.begin(), std::next(ranked.begin(), best)) == recorded);
    constexpr auto row_7 = std::size_t{6};
    CHECK(recorded.at(row_7) == Positions{260, 450});
}
