#include "core/analysis/legacy_verdict.hpp"
#include "core/io/genbank.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <charconv>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace {

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

auto number(std::string const& text) -> std::size_t {
    auto value = std::size_t{0};
    auto const* const last = std::next(text.data(), std::ssize(text));
    auto const [end, error] = std::from_chars(text.data(), last, value);
    if (error != std::errc{} or end != last) {
        throw std::runtime_error{"not a number: " + text};
    }
    return value;
}

// columns: settings, molecule, target_a, start_offset_a, target_b,
// end_offset_b, length, outcome
auto recorded_searches() -> std::vector<FindPairSettings> {
    constexpr auto column_count = 8U;
    auto searches = std::vector<FindPairSettings>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "search" / "find_pair_no_match.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        if (fields.size() != column_count or fields[0] != "built-in" or fields[1] != "TLONG.GB" or
            fields[7] != "no_matching_primers") {
            throw std::runtime_error{"unexpected line: " + line};
        }
        searches.push_back({.target_a = number(fields[2]),
                            .start_offset_a = number(fields[3]),
                            .target_b = number(fields[4]),
                            .end_offset_b = number(fields[5]),
                            .length = number(fields[6])});
    }
    return searches;
}

}  // namespace

// Reference values: PRIMER.exe, "No Matching Primers Found" instead of
// the Find Primer Pair list (tests/legacy/search/find_pair_no_match.tsv,
// recorded 2026-10-19)

TEST_CASE("Find Pair: every recorded No Matching Primers Found has no pair", "[search][legacy]") {
    auto const molecule = tlong();
    auto const searches = recorded_searches();
    REQUIRE(searches.size() == 2);
    for (auto const& search : searches) {
        CAPTURE(search.target_a, search.target_b);
        auto const result = legacy_find_pair(molecule, search, primerpicker::analysis::LegacyCriteria{});
        CHECK(result.pairs.empty());
        CHECK(result.limit == FindPairLimit::none);
    }
}
