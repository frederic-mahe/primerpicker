#include "core/analysis/legacy_verdict.hpp"
#include "core/io/genbank.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <iterator>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using primerpicker::analysis::LegacyActualValues;
using primerpicker::analysis::LegacyCriteria;
using primerpicker::search::FindPairSettings;
using primerpicker::search::FindPairWeights;
using primerpicker::search::legacy_find_pair;
using primerpicker::search::legacy_pair_score;
using primerpicker::search::legacy_rank_pairs;
using primerpicker::search::PairCandidate;
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

auto legacy_path(std::string const& directory, std::string const& name) -> std::filesystem::path {
    return std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / directory / name;
}

auto molecule(std::filesystem::path const& path) -> std::vector<Nucleotide> {
    auto const parsed = primerpicker::io::parse_genbank(read_file(path));
    auto const sequence = parsed.molecule.has_value()
                              ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                              : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    return *sequence;
}

auto split(std::string const& line, char const separator) -> std::vector<std::string> {
    auto fields = std::vector<std::string>{};
    auto stream = std::istringstream{line};
    for (auto field = std::string{}; std::getline(stream, field, separator);) {
        fields.push_back(field);
    }
    return fields;
}

// the data lines of a recorded table, as fields
auto table(std::filesystem::path const& path) -> std::vector<std::vector<std::string>> {
    auto rows = std::vector<std::vector<std::string>>{};
    auto file = std::ifstream{path};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (not line.starts_with('#')) {
            rows.push_back(split(line, '\t'));
        }
    }
    return rows;
}

// "0 5 0 0 1": GC mismatch, Tm mismatch, 3' dimers, any dimers, false priming
auto weights(std::string const& text) -> FindPairWeights {
    auto const digits = split(text, ' ');
    if (digits.size() != 5) {
        throw std::runtime_error{"five weights expected: " + text};
    }
    auto const weight = [&digits](std::size_t const index) { return std::stoi(digits.at(index)); };
    return FindPairWeights{.gc_mismatch = weight(0),
                           .tm_mismatch = weight(1),
                           .dimer_3prime = weight(2),
                           .dimer_any = weight(3),
                           .false_priming = weight(4)};
}

auto tm_text(std::optional<int> const site_tm) -> std::string {
    return site_tm.has_value() ? std::to_string(*site_tm) : std::string{"-"};
}

// the values shown on a line of the list, as the recorded tables write them
auto describe(PairCandidate const& pair) -> std::string {
    return std::format("{} {} {} {} {} {} {} {} {} {}", pair.position_a, pair.position_b, pair.primer_a.gc_percent,
                       pair.primer_b.gc_percent, pair.primer_a.tm_celsius, pair.primer_b.tm_celsius,
                       pair.max_dimer_3prime_matches, pair.max_dimer_adjacent_bases,
                       tm_text(pair.false_priming_tm_a), tm_text(pair.false_priming_tm_b));
}

auto positions(PairCandidate const& pair) -> std::string {
    return std::format("{} {}", pair.position_a, pair.position_b);
}

auto ranked(std::vector<Nucleotide> const& sequence, FindPairSettings const& settings, FindPairWeights const& weights)
    -> std::vector<PairCandidate> {
    return legacy_rank_pairs(legacy_find_pair(sequence, settings, LegacyCriteria{}).pairs, weights);
}

auto join(std::vector<std::string> const& fields, std::size_t const first) -> std::string {
    auto text = std::string{};
    for (auto index = first; index < fields.size(); ++index) {
        text += (text.empty() ? "" : " ") + fields[index];
    }
    return text;
}

constexpr auto no_site = std::optional<int>{};

// the values of a pair that enter its score (dimers: 2 and 4)
struct ScoreTerms {
    int gc_difference = 0;
    int tm_difference = 0;
    std::optional<int> false_priming_tm_a;
    std::optional<int> false_priming_tm_b;
};

auto candidate(ScoreTerms const& terms) -> PairCandidate {
    constexpr auto gc_percent = 50;
    constexpr auto tm_celsius = 60;
    auto primer_a = LegacyActualValues{};
    primer_a.gc_percent = gc_percent;
    primer_a.tm_celsius = tm_celsius;
    auto primer_b = primer_a;
    primer_b.gc_percent += terms.gc_difference;
    primer_b.tm_celsius -= terms.tm_difference;
    return PairCandidate{.position_a = 1,
                         .position_b = 2,
                         .primer_a = primer_a,
                         .primer_b = primer_b,
                         .max_dimer_3prime_matches = 2,
                         .max_dimer_adjacent_bases = 4,
                         .false_priming_tm_a = terms.false_priming_tm_a,
                         .false_priming_tm_b = terms.false_priming_tm_b};
}

}  // namespace

// Reference values: PRIMER.exe, Find Primer Pair lists
// (tests/legacy/search/find_pair.tsv, find_pair_weights.tsv and
// find_pair_false_priming.tsv, recorded 2026-10-18; README, "Find Pair")

TEST_CASE("Find Pair ranks the pairs as the recorded list (built-in weights)", "[search][legacy]") {
    auto const computed = ranked(molecule(legacy_path("pairs", "TEST600.GB")),
                                 FindPairSettings{.target_a = 150, .target_b = 400}, FindPairWeights{});
    auto expected = std::vector<std::string>{};
    constexpr auto first_value_column = 7U;
    std::ranges::transform(table(legacy_path("search", "find_pair.tsv")), std::back_inserter(expected),
                           [](std::vector<std::string> const& fields) { return join(fields, first_value_column); });
    auto actual = std::vector<std::string>{};
    std::ranges::transform(computed, std::back_inserter(actual), describe);
    REQUIRE(expected.size() == 10);
    CHECK(actual == expected);
}

TEST_CASE("Find Pair ranks the pairs as the recorded lists (other weights)", "[search][legacy]") {
    auto const sequence = molecule(legacy_path("pairs", "TEST600.GB"));
    auto expected = std::map<std::string, std::vector<std::string>>{};  // weights -> positions, by rank
    for (auto const& fields : table(legacy_path("search", "find_pair_weights.tsv"))) {
        expected[fields.at(0)].push_back(fields.at(2) + " " + fields.at(3));
    }
    REQUIRE(expected.size() == 2);
    for (auto const& [weight_text, expected_positions] : expected) {
        auto const computed =
            ranked(sequence, FindPairSettings{.target_a = 150, .target_b = 400}, weights(weight_text));
        auto actual = std::vector<std::string>{};
        std::ranges::transform(computed, std::back_inserter(actual), positions);
        INFO("weights " << weight_text);
        CHECK(actual == expected_positions);
    }
}

TEST_CASE("Find Pair ranks pairs with false-priming sites as the recorded lists", "[search][legacy]") {
    // (molecule, weights) -> pair values, by rank
    auto expected = std::map<std::pair<std::string, std::string>, std::vector<std::string>>{};
    constexpr auto first_value_column = 5U;
    for (auto const& fields : table(legacy_path("search", "find_pair_false_priming.tsv"))) {
        expected[{fields.at(0), fields.at(3)}].push_back(join(fields, first_value_column));
    }
    REQUIRE(expected.size() == 9);  // TFP.GB: 7 weight sets, TFP2.GB: 2
    for (auto const& [run, expected_pairs] : expected) {
        auto const& [molecule_name, weight_text] = run;
        auto const computed = ranked(molecule(legacy_path("search", molecule_name)),
                                     FindPairSettings{.target_a = 150, .target_b = 400}, weights(weight_text));
        auto actual = std::vector<std::string>{};
        std::ranges::transform(computed, std::back_inserter(actual), describe);
        INFO(molecule_name << ", weights " << weight_text);
        REQUIRE(expected_pairs.size() == 30);  // PRIMER.exe: "Pairs: 30 possible"
        CHECK(actual == expected_pairs);
    }
}

// The cases below restate, on single pairs, the score that the lists
// above follow (tests/legacy/README.md, "Find Pair")

TEST_CASE("Find Pair score: false priming counts the site Tm of both primers", "[search][legacy]") {
    constexpr auto false_priming_only = FindPairWeights{
        .gc_mismatch = 0, .tm_mismatch = 0, .dimer_3prime = 0, .dimer_any = 0, .false_priming = 1};
    // TFP.GB: 124/471, whose only site (primer B) has a Tm of 6, scores 6
    auto const site_b = candidate(
        {.gc_difference = 0, .tm_difference = 0, .false_priming_tm_a = no_site, .false_priming_tm_b = 6});
    CHECK(legacy_pair_score(site_b, false_priming_only) == 6);
    // TFP.GB, weights 0 1 0 0 1: 78/471 (Tm difference 2, sites of Tm 32
    // and 6) ranks between 79/460 (score 0 + 39) and 79/461 (2 + 39)
    auto const sites_a_and_b =
        candidate({.gc_difference = 0, .tm_difference = 0, .false_priming_tm_a = 32, .false_priming_tm_b = 6});
    CHECK(legacy_pair_score(sites_a_and_b, false_priming_only) == 38);
    auto const no_sites = candidate(
        {.gc_difference = 0, .tm_difference = 0, .false_priming_tm_a = no_site, .false_priming_tm_b = no_site});
    CHECK(legacy_pair_score(no_sites, false_priming_only) == 0);
}

TEST_CASE("Find Pair score: weighted sum of the differences and dimers", "[search][legacy]") {
    // TEST600.GB, built-in weights: 73/425 (GC 50/50, Tm 69/67, dimers
    // 2 and 5) comes first (find_pair.tsv); here, dimers 2 and 4
    auto const tm_difference_2 = candidate(
        {.gc_difference = 0, .tm_difference = 2, .false_priming_tm_a = no_site, .false_priming_tm_b = no_site});
    CHECK(legacy_pair_score(tm_difference_2, FindPairWeights{}) == 0 + 2 + 2 + 4);
    constexpr auto tm_times_five = FindPairWeights{
        .gc_mismatch = 0, .tm_mismatch = 5, .dimer_3prime = 0, .dimer_any = 0, .false_priming = 1};
    // TEST600.GB, weights 0 5 0 0 1: 71/453 (Tm 75/69, a site of Tm 21)
    // comes after pairs of score 45 (find_pair_weights.tsv)
    auto const site_a = candidate(
        {.gc_difference = 5, .tm_difference = 6, .false_priming_tm_a = 21, .false_priming_tm_b = no_site});
    CHECK(legacy_pair_score(site_a, tm_times_five) == 30 + 21);
}

TEST_CASE("Find Pair ranking keeps generation order between equal scores", "[search][legacy]") {
    // TEST600.GB, weights 0 0 0 0 1: the eight pairs without
    // false-priming site keep their generation order (find_pair_weights.tsv)
    auto pairs = std::vector<PairCandidate>{
        candidate({.gc_difference = 0, .tm_difference = 0, .false_priming_tm_a = 21, .false_priming_tm_b = no_site}),
        candidate(
            {.gc_difference = 0, .tm_difference = 0, .false_priming_tm_a = no_site, .false_priming_tm_b = no_site}),
        candidate(
            {.gc_difference = 5, .tm_difference = 9, .false_priming_tm_a = no_site, .false_priming_tm_b = no_site})};
    pairs[0].position_a = 71;
    pairs[1].position_a = 72;
    pairs[2].position_a = 73;
    constexpr auto false_priming_only = FindPairWeights{
        .gc_mismatch = 0, .tm_mismatch = 0, .dimer_3prime = 0, .dimer_any = 0, .false_priming = 1};
    auto const result = legacy_rank_pairs(std::move(pairs), false_priming_only);
    REQUIRE(result.size() == 3);
    CHECK(result[0].position_a == 72);
    CHECK(result[1].position_a == 73);
    CHECK(result[2].position_a == 71);
}
