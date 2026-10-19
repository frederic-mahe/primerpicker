#include "core/analysis/legacy_verdict.hpp"
#include "core/io/genbank.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/thermo/units.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <ios>
#include <iterator>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using primerpicker::analysis::LegacyCriteria;
using primerpicker::search::FindPairLimit;
using primerpicker::search::FindPairSettings;
using primerpicker::search::legacy_find_pair;
using primerpicker::sequence::Nucleotide;

auto read_file(std::filesystem::path const& path) -> std::string {
    auto text = std::string(std::filesystem::file_size(path), '\0');
    auto file = std::ifstream{path, std::ios::binary};
    file.read(text.data(), std::ssize(text));
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    return text;
}

auto tlong() -> std::vector<Nucleotide> {
    auto const parsed = primerpicker::io::parse_genbank(
        read_file(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "enzymes" / "TLONG.GB"));
    auto const sequence = parsed.molecule.has_value()
                              ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                              : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot read TLONG.GB"};
    }
    return *sequence;
}

// tests/legacy/config/PD2_RELAXED.CFG, read with the layout of
// tests/legacy/README.md ("Settings file pd2.cfg"): GC 0-100 %, Tm
// 0-99 °C, hairpin cutoff -99.0 kcal, every reject threshold 99
constexpr auto relaxed_criteria = LegacyCriteria{
    .gc_min_percent = 0,
    .gc_max_percent = 100,
    .tm_min_celsius = 0,
    .tm_max_celsius = 99,
    .hairpin_cutoff = primerpicker::thermo::KilocaloriesPerMole{-99.0},
    .reject_dimer_3prime_matches = 99,
    .reject_dimer_adjacent_bases = 99,
    .reject_run_length = 99,
    .reject_gc_run_3prime = 99,
};

auto criteria(std::string_view const settings) -> LegacyCriteria {
    if (settings == "PD2_RELAXED.CFG") {
        return relaxed_criteria;
    }
    if (settings == "built-in") {
        return LegacyCriteria{};
    }
    throw std::runtime_error{"unknown settings: " + std::string{settings}};
}

auto number(std::string const& text) -> std::size_t {
    auto value = std::size_t{0};
    auto const* const last = std::next(text.data(), std::ssize(text));
    auto const [end, error] = std::from_chars(text.data(), last, value);
    if (error != std::errc{} or end != last) {
        throw std::runtime_error{"not a number: " + text};
    }
    return value;
}

// "-": not shown by PRIMER.exe
auto shown(std::string const& text) -> std::optional<std::size_t> {
    return text == "-" ? std::nullopt : std::optional{number(text)};
}

struct Recording {
    std::string settings;
    FindPairSettings search;
    std::optional<std::size_t> possible_a;
    std::optional<std::size_t> possible_b;
    std::optional<std::size_t> pairs;
    std::string outcome;
};

// columns: settings, molecule, target_a, start_offset_a, target_b,
// end_offset_b, length, possible_a, possible_b, pairs, outcome
auto recordings() -> std::vector<Recording> {
    constexpr auto column_count = 11U;
    auto rows = std::vector<Recording>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "search" / "find_pair_limits.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        if (fields.size() != column_count or fields[1] != "TLONG.GB") {
            throw std::runtime_error{"unexpected line: " + line};
        }
        rows.push_back({.settings = fields[0],
                        .search = {.target_a = number(fields[2]),
                                   .start_offset_a = number(fields[3]),
                                   .target_b = number(fields[4]),
                                   .end_offset_b = number(fields[5]),
                                   .length = number(fields[6])},
                        .possible_a = shown(fields[7]),
                        .possible_b = shown(fields[8]),
                        .pairs = shown(fields[9]),
                        .outcome = fields[10]});
    }
    return rows;
}

auto expected_limit(std::string_view const outcome) -> FindPairLimit {
    if (outcome == "too_many_sites") {
        return FindPairLimit::too_many_primers;
    }
    if (outcome == "pairs_capped") {
        return FindPairLimit::pairs_capped;
    }
    return FindPairLimit::none;
}

// the counts PRIMER.exe showed, and the outcome, of one recorded search
void check_recording(std::vector<Nucleotide> const& molecule, Recording const& row) {
    CAPTURE(row.settings, row.search.target_a, row.search.start_offset_a, row.search.target_b,
            row.search.end_offset_b, row.outcome);
    auto const result = legacy_find_pair(molecule, row.search, criteria(row.settings));
    CHECK(result.limit == expected_limit(row.outcome));
    CHECK(row.possible_a.value_or(result.possible_a) == result.possible_a);
    CHECK(row.possible_b.value_or(result.possible_b) == result.possible_b);
    CHECK(row.pairs.value_or(result.pairs.size()) == result.pairs.size());
    auto const no_pair_expected = row.outcome == "too_many_sites" or row.outcome == "no_matching_primers";
    CHECK((result.pairs.empty() or not no_pair_expected));
}

}  // namespace

// Reference values: PRIMER.exe, Find Primer Pair on TLONG.GB
// (tests/legacy/search/find_pair_limits.tsv, recorded 2026-10-19)

TEST_CASE("Find Pair limits: every recorded search", "[search][legacy]") {
    auto const molecule = tlong();
    auto const rows = recordings();
    REQUIRE(rows.size() == 11);
    for (auto const& row : rows) {
        check_recording(molecule, row);
    }
}

TEST_CASE("Find Pair limits: 99 possible primers are searched", "[search][legacy]") {
    // PRIMER.exe (find_pair_limits.tsv): relaxed settings, target A 400,
    // offset 98: "Primer A: 99 possible primers", "Primer B: 6", "Pairs: 220"
    auto const result = legacy_find_pair(
        tlong(), {.target_a = 400, .start_offset_a = 98, .target_b = 800, .end_offset_b = 5}, relaxed_criteria);
    CHECK(result.limit == FindPairLimit::none);
    CHECK(result.possible_a == 99);
    CHECK(result.pairs.size() == 220);
}

TEST_CASE("Find Pair limits: 100 possible primers stop the search", "[search][legacy]") {
    // PRIMER.exe (find_pair_limits.tsv): relaxed settings, target A 400,
    // offset 99: "More than 100 sites found.  Suggest you redefine
    // matching limits.", and no search
    auto const result = legacy_find_pair(
        tlong(), {.target_a = 400, .start_offset_a = 99, .target_b = 800, .end_offset_b = 5}, relaxed_criteria);
    CHECK(result.limit == FindPairLimit::too_many_primers);
    CHECK(result.pairs.empty());
}

TEST_CASE("Find Pair limits: the limit counts possible primers, not positions", "[search][legacy]") {
    // PRIMER.exe (find_pair_limits.tsv): built-in criteria, 111
    // positions for primer A: "Primer A: 9 possible primers", "Primer B:
    // 16 possible primers", "Pairs: 48 possible"
    auto const result = legacy_find_pair(
        tlong(), {.target_a = 200, .start_offset_a = 110, .target_b = 800, .end_offset_b = 50}, LegacyCriteria{});
    CHECK(result.limit == FindPairLimit::none);
    CHECK(result.possible_a == 9);
    CHECK(result.possible_b == 16);
    CHECK(result.pairs.size() == 48);
}

TEST_CASE("Find Pair limits: the pairs stop at 1000, in generation order", "[search][legacy]") {
    // PRIMER.exe (find_pair_limits.tsv): relaxed settings, 61 x 61
    // positions: the count stops at "Pairs:  1000 possible" (then the
    // program hangs; version 1 keeps those 1000 pairs and warns)
    auto const result = legacy_find_pair(
        tlong(), {.target_a = 400, .start_offset_a = 60, .target_b = 800, .end_offset_b = 60}, relaxed_criteria);
    CHECK(result.limit == FindPairLimit::pairs_capped);
    CHECK(result.possible_a == 61);
    CHECK(result.possible_b == 61);
    REQUIRE(result.pairs.size() == 1000);
    auto const by_generation = [](auto const& pair) { return std::pair{pair.position_a, pair.position_b}; };
    CHECK(std::ranges::is_sorted(result.pairs, std::less{}, by_generation));
}
