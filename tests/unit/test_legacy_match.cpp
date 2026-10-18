#include "core/analysis/legacy_verdict.hpp"
#include "core/io/genbank.hpp"
#include "core/search/legacy_match.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <iterator>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using primerpicker::search::legacy_match_active_primer;
using primerpicker::search::MatchCandidate;
using primerpicker::search::MatchSettings;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::primer_at;
using primerpicker::sequence::Strand;

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

// the columns of match_active.tsv, from position on
auto describe(MatchCandidate const& candidate) -> std::string {
    auto const& actual = candidate.actual;
    return std::format("{} {} {} {} {} {} {} {} {} {}", candidate.position, candidate.distance, actual.gc_percent,
                       actual.tm_celsius, actual.dimer_3prime_matches, actual.dimer_adjacent_bases,
                       candidate.inter_dimer_3prime_matches, candidate.inter_dimer_adjacent_bases, actual.longest_run,
                       actual.gc_run_3prime);
}

}  // namespace

// Reference values: PRIMER.exe, Find Matching Primer list
// (tests/legacy/search/match_active.tsv, recorded 2026-10-18)

TEST_CASE("primers are taken from either strand of a molecule", "[sequence]") {
    auto const molecule = test600();
    // PRIMER.exe: position 51 -> CTCAAAGGCTCAGACATCAG, 400(C) -> TGTATCCAGGTAATGCAAGG
    CHECK(primer_at(molecule, {.position = 51, .length = 20, .strand = Strand::top}) ==
          primerpicker::sequence::to_nucleotides("CTCAAAGGCTCAGACATCAG"));
    CHECK(primer_at(molecule, {.position = 400, .length = 20, .strand = Strand::complement}) ==
          primerpicker::sequence::to_nucleotides("TGTATCCAGGTAATGCAAGG"));
    // beyond the molecule
    CHECK_FALSE(primer_at(molecule, {.position = 590, .length = 20, .strand = Strand::top}).has_value());
    CHECK_FALSE(primer_at(molecule, {.position = 10, .length = 20, .strand = Strand::complement}).has_value());
}

TEST_CASE("Match Active Primer finds the recorded candidates", "[search][legacy]") {
    auto const molecule = test600();
    auto const primer_a = primer_at(molecule, {.position = 51, .length = 20, .strand = Strand::top});
    REQUIRE(primer_a.has_value());
    auto const candidates =
        legacy_match_active_primer({.molecule = molecule, .primer_a = primer_a.value_or(std::vector<Nucleotide>{})},
                                   MatchSettings{.target = 400}, primerpicker::analysis::LegacyCriteria{});
    auto computed = std::vector<std::string>{};
    std::ranges::transform(candidates, std::back_inserter(computed), describe);

    auto expected = std::vector<std::string>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "search" / "match_active.tsv"};
    constexpr auto first_value_column = 8U;
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
        for (auto const& field : std::span{fields}.subspan(first_value_column)) {
            values += (values.empty() ? "" : " ") + field;
        }
        expected.push_back(values);
    }
    REQUIRE(expected.size() == 11);
    CHECK(computed == expected);
}
