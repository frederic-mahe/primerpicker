#include "core/io/molecule_file.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using primerpicker::io::LegacySystemDefaults;
using primerpicker::sequence::Nucleotide;


auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

auto read_file(std::filesystem::path const& path) -> std::string {
    auto file = std::ifstream{path, std::ios::binary};
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    auto text = std::ostringstream{};
    text << file.rdbuf();
    return std::move(text).str();
}

auto defaults_from(std::string_view const config) -> LegacySystemDefaults {
    if (config == "-") {
        return LegacySystemDefaults{};
    }
    auto const text = read_file(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "config" / config);
    auto bytes = std::vector<std::byte>{};
    std::ranges::transform(text, std::back_inserter(bytes), [](char const character) {
        return static_cast<std::byte>(character);
    });
    auto const result = primerpicker::io::parse_pd2_cfg(bytes);
    if (not result.defaults.has_value()) {
        throw std::runtime_error{"cannot parse " + std::string{config}};
    }
    return *result.defaults;
}

auto number(std::string const& text) -> int {
    auto value = 0;
    auto const* const last = std::next(text.data(), std::ssize(text));
    auto const [end, error] = std::from_chars(text.data(), last, value);
    if (error != std::errc{} or end != last) {
        throw std::runtime_error{"not a number: " + text};
    }
    return value;
}

// search/find_pair_3prime_runs.tsv: the possible primers and pairs
// counted by PRIMER.exe on search/TRUN.GB (offsets 100, length 20)
struct RecordedCounts {
    std::string config;
    std::size_t target_a;
    std::size_t target_b;
    std::size_t possible_a;
    std::size_t possible_b;
    std::size_t pairs;
    std::string note;
};

auto recorded_counts() -> std::vector<RecordedCounts> {
    constexpr auto column_count = 7U;
    auto counts = std::vector<RecordedCounts>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "search" / "find_pair_3prime_runs.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        if (fields.size() != column_count) {
            throw std::runtime_error{"unexpected line: " + line};
        }
        auto const count = [](std::string const& text) { return static_cast<std::size_t>(number(text)); };
        counts.push_back({.config = fields[0],
                          .target_a = count(fields[1]),
                          .target_b = count(fields[2]),
                          .possible_a = count(fields[3]),
                          .possible_b = count(fields[4]),
                          .pairs = count(fields[5]),
                          .note = fields[6]});
    }
    return counts;
}

auto trun_molecule() -> std::vector<Nucleotide> {
    auto const text = read_file(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "search" / "TRUN.GB");
    auto const parsed = primerpicker::io::parse_molecule_file({.text = text, .name = "TRUN.GB"});
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"cannot parse TRUN.GB"};
    }
    return nucleotides(parsed.molecule->sequence);
}

auto search(RecordedCounts const& recorded, std::vector<Nucleotide> const& molecule)
    -> primerpicker::search::FindPairResult {
    auto const defaults = defaults_from(recorded.config);
    auto const settings = primerpicker::search::FindPairSettings{.target_a = recorded.target_a,
                                                                 .start_offset_a = defaults.find_pair_start_offset,
                                                                 .target_b = recorded.target_b,
                                                                 .end_offset_b = defaults.find_pair_end_offset,
                                                                 .length = defaults.primer_length,
                                                                 .gc_tolerance_percent = defaults.gc_tolerance_percent,
                                                                 .tm_tolerance_celsius = defaults.tm_tolerance_celsius};
    return primerpicker::search::legacy_find_pair(molecule, settings, defaults.criteria);
}

}  // namespace

// Sweep (tests/legacy/README.md, "Find Pair and runs at the 3' end"):
// Find Pair rejects a primer whose 3'-end run reaches the runs criterion,
// a run that the Primer Analysis screen does not count
TEST_CASE("Find Pair: possible primers, 3'-end runs included", "[search][legacy]") {
    auto const molecule = trun_molecule();
    auto const counts = recorded_counts();
    REQUIRE(counts.size() >= 8);
    for (auto const& recorded : counts) {
        INFO(recorded.config << " " << recorded.target_a << " " << recorded.target_b);
        auto const result = search(recorded, molecule);
        CHECK(result.possible_a == recorded.possible_a);
        if (recorded.note == "full-length-copy") {
            continue;  // primers B with a copy elsewhere in the molecule (README): A only
        }
        CHECK(result.possible_b == recorded.possible_b);
        CHECK(result.pairs.size() == recorded.pairs);
    }
}
