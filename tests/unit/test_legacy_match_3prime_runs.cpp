#include "core/io/molecule_file.hpp"
#include "core/search/legacy_match.hpp"
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

auto number(std::string const& text) -> int {
    auto value = 0;
    auto const* const last = std::next(text.data(), std::ssize(text));
    auto const [end, error] = std::from_chars(text.data(), last, value);
    if (error != std::errc{} or end != last) {
        throw std::runtime_error{"not a number: " + text};
    }
    return value;
}

// search/match_3prime_runs.tsv: the Find Matching Primer lists printed
// by PRIMER.exe (reports/match/), form defaults
struct RecordedList {
    std::string molecule;  // relative to tests/legacy
    std::size_t primer_a_position;
    std::size_t target;
    std::vector<std::size_t> positions;
};

auto recorded_lists() -> std::vector<RecordedList> {
    constexpr auto column_count = 4U;
    auto lists = std::vector<RecordedList>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "search" / "match_3prime_runs.tsv"};
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
        auto positions = std::vector<std::size_t>{};
        auto position_stream = std::istringstream{fields[3]};
        for (auto position = std::string{}; std::getline(position_stream, position, ',');) {
            positions.push_back(static_cast<std::size_t>(number(position)));
        }
        lists.push_back({.molecule = fields[0],
                         .primer_a_position = static_cast<std::size_t>(number(fields[1])),
                         .target = static_cast<std::size_t>(number(fields[2])),
                         .positions = positions});
    }
    return lists;
}

auto molecule_from(std::string const& path) -> std::vector<Nucleotide> {
    auto const text = read_file(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / path);
    auto const parsed = primerpicker::io::parse_molecule_file({.text = text, .name = path});
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"cannot parse " + path};
    }
    return nucleotides(parsed.molecule->sequence);
}

}  // namespace

// Sweep (tests/legacy/README.md, "Find Pair and runs at the 3' end"):
// Match Active Primer, like Find Pair, rejects a primer B whose 3'-end
// run reaches the runs criterion
TEST_CASE("Match Active Primer: the recorded lists, 3'-end runs included", "[search][legacy]") {
    auto const lists = recorded_lists();
    REQUIRE(lists.size() >= 2);
    for (auto const& recorded : lists) {
        INFO(recorded.molecule << " " << recorded.primer_a_position << " " << recorded.target);
        auto const molecule = molecule_from(recorded.molecule);
        auto const primer_a = std::span{molecule}.subspan(recorded.primer_a_position - 1, 20);
        auto const list = primerpicker::search::legacy_match_active_primer(
            {.molecule = molecule, .primer_a = primer_a}, {.target = recorded.target},
            primerpicker::analysis::LegacyCriteria{});
        auto positions = std::vector<std::size_t>{};
        std::ranges::transform(list, std::back_inserter(positions),
                               [](primerpicker::search::MatchCandidate const& candidate) { return candidate.position; });
        CHECK(positions == recorded.positions);
    }
}
