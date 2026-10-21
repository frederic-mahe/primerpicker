#include "core/io/molecule_file.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <filesystem>
#include <fstream>
#include <ios>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using primerpicker::search::FindPairLimit;
using primerpicker::search::FindPairSettings;
using primerpicker::search::legacy_find_pair;
using primerpicker::sequence::Nucleotide;

auto test600() -> std::vector<Nucleotide> {
    auto const path = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "pairs" / "TEST600.GB";
    auto file = std::ifstream{path, std::ios::binary};
    auto text = std::ostringstream{};
    text << file.rdbuf();
    auto const parsed = primerpicker::io::parse_molecule_file({.text = text.str(), .name = "TEST600.GB"});
    auto sequence = parsed.molecule.has_value() ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                                                : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot parse TEST600.GB"};
    }
    return std::move(*sequence);
}

constexpr auto settings = FindPairSettings{.target_a = 150, .target_b = 400};

}  // namespace

// PRIMER.exe's "Press key to interrupt search" (from the strings of
// examples/PRIMER.exe): an interrupted search returns no pair
TEST_CASE("Find Pair: an interrupted search returns no pair", "[search]") {
    auto const interrupt = std::atomic<bool>{true};
    auto const result = legacy_find_pair(test600(), settings, {}, {}, interrupt);
    CHECK(result.limit == FindPairLimit::interrupted);
    CHECK(result.pairs.empty());
}

TEST_CASE("Find Pair: a search that is not interrupted is the usual one", "[search]") {
    auto const interrupt = std::atomic<bool>{false};
    auto const molecule = test600();
    auto const result = legacy_find_pair(molecule, settings, {}, {}, interrupt);
    auto const usual = legacy_find_pair(molecule, settings, {});
    CHECK(result.limit == FindPairLimit::none);
    CHECK(result.pairs.size() == usual.pairs.size());
    CHECK(result.possible_a == usual.possible_a);
    CHECK(result.possible_b == usual.possible_b);
}
