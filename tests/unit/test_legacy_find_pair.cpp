#include "core/analysis/legacy_verdict.hpp"
#include "core/io/genbank.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <iterator>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using primerpicker::search::FindPairSettings;
using primerpicker::search::legacy_find_pair;
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

auto test600() -> std::vector<Nucleotide> {
    auto const parsed =
        primerpicker::io::parse_genbank(read_file(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "pairs" / "TEST600.GB"));
    auto const sequence = parsed.molecule.has_value()
                              ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                              : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot read TEST600.GB"};
    }
    return *sequence;
}

auto tm_text(std::optional<int> const site_tm) -> std::string {
    return site_tm.has_value() ? std::to_string(*site_tm) : std::string{"-"};
}

// the columns of find_pair.tsv, from position_a on
auto describe(PairCandidate const& pair) -> std::string {
    return std::format("{} {} {} {} {} {} {} {} {} {}", pair.position_a, pair.position_b, pair.primer_a.gc_percent,
                       pair.primer_b.gc_percent, pair.primer_a.tm_celsius, pair.primer_b.tm_celsius,
                       pair.max_dimer_3prime_matches, pair.max_dimer_adjacent_bases,
                       tm_text(pair.false_priming_tm_a), tm_text(pair.false_priming_tm_b));
}

// the recorded pairs, as describe() writes them
auto recorded_pairs() -> std::vector<std::string> {
    auto pairs = std::vector<std::string>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "search" / "find_pair.tsv"};
    constexpr auto first_value_column = 7U;
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        auto values = std::string{};
        for (auto index = first_value_column; index < fields.size(); ++index) {
            values += (values.empty() ? "" : " ") + fields[index];
        }
        pairs.push_back(values);
    }
    return pairs;
}

}  // namespace

// Reference values: PRIMER.exe, Find Primer Pair list
// (tests/legacy/search/find_pair.tsv, recorded 2026-10-18)

TEST_CASE("Find Pair finds the recorded pairs and their values", "[search][legacy]") {
    auto const molecule = test600();
    auto const result = legacy_find_pair(molecule, FindPairSettings{.target_a = 150, .target_b = 400},
                                         primerpicker::analysis::LegacyCriteria{});
    CHECK(result.possible_a == 3);  // PRIMER.exe: "Primer A: 3 possible primers"
    CHECK(result.possible_b == 7);  // PRIMER.exe: "Primer B: 7 possible primers"
    auto computed = std::vector<std::string>{};
    std::ranges::transform(result.pairs, std::back_inserter(computed), describe);
    std::ranges::sort(computed);

    auto expected = recorded_pairs();
    std::ranges::sort(expected);
    REQUIRE(expected.size() == 10);  // PRIMER.exe: "Pairs: 10 possible"
    CHECK(computed == expected);
}
