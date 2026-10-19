#include "core/analysis/legacy_verdict.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/search/legacy_match.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using primerpicker::analysis::LegacyCriteria;
using primerpicker::search::FindPairSettings;
using primerpicker::search::legacy_find_pair;
using primerpicker::search::legacy_match_active_primer;
using primerpicker::search::MatchSettings;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::primer_at;
using primerpicker::sequence::Strand;
using primerpicker::sequence::window_end;

constexpr auto largest = std::numeric_limits<std::size_t>::max();

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return std::move(*converted);
}

// 60 bases, primers of 20 fit anywhere on it
auto molecule() -> std::vector<Nucleotide> {
    return nucleotides("GTAAAACGACGGCCAGTGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGG");
}

}  // namespace

// Positions far beyond the molecule (found by fuzz_cli, 2026-10-19):
// an unsigned overflow made primer_at read outside the molecule, and a
// search window ending at the largest size_t made a search loop forever.
// A primer that does not fit on the molecule does not exist.

TEST_CASE("primer_at: positions beyond the molecule", "[sequence][limits]") {
    auto const bases = molecule();
    CHECK(not primer_at(bases, {.position = largest, .length = 20, .strand = Strand::top}).has_value());
    CHECK(not primer_at(bases, {.position = largest - 10, .length = 20, .strand = Strand::top}).has_value());
    CHECK(not primer_at(bases, {.position = 61, .length = 1, .strand = Strand::top}).has_value());
    CHECK(not primer_at(bases, {.position = 42, .length = largest, .strand = Strand::top}).has_value());
    CHECK(not primer_at(bases, {.position = largest, .length = 20, .strand = Strand::complement}).has_value());
    CHECK(primer_at(bases, {.position = 41, .length = 20, .strand = Strand::top}).has_value());  // the last fit
}

TEST_CASE("window_end: the end of a search window, within the molecule", "[sequence][limits]") {
    CHECK(window_end(400, 100, 600) == 500);
    CHECK(window_end(550, 100, 600) == 600);
    CHECK(window_end(largest - 100, 100, 600) == 600);
    CHECK(window_end(largest, largest, 600) == 600);
    CHECK(window_end(1, 0, 0) == 0);
}

TEST_CASE("Find Pair: targets beyond the molecule find no primer", "[search][limits]") {
    auto const bases = molecule();
    auto const result =
        legacy_find_pair(bases, FindPairSettings{.target_a = 16, .target_b = largest - 100}, LegacyCriteria{});
    CHECK(result.possible_b == 0);
    CHECK(result.pairs.empty());
    auto const far_a = legacy_find_pair(bases, FindPairSettings{.target_a = largest, .target_b = 40}, LegacyCriteria{});
    CHECK(far_a.possible_a == 0);
}

TEST_CASE("Match Active Primer: a target beyond the molecule finds no primer", "[search][limits]") {
    auto const bases = molecule();
    auto const primer_a = nucleotides("GTAAAACGACGGCCAGT");
    auto const candidates = legacy_match_active_primer({.molecule = bases, .primer_a = primer_a},
                                                       MatchSettings{.target = largest - 50}, LegacyCriteria{});
    CHECK(candidates.empty());
}
