#include "core/enzymes/site_search.hpp"
#include "core/io/molecule_file.hpp"
#include "core/search/legacy_find_sequence.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/topology.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using primerpicker::search::legacy_find_sequence;
using primerpicker::search::legacy_site_list_lines;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::Topology;

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

auto test600() -> std::vector<Nucleotide> {
    auto const parsed = primerpicker::io::parse_molecule_file(
        {.text = read_file(legacy_dir / "pairs" / "TEST600.GB"), .name = "TEST600.GB"});
    auto sequence = parsed.molecule.has_value() ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                                                : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot parse TEST600.GB"};
    }
    return std::move(*sequence);
}

// The site list of the last FIND SEQUENCE screen of a recording: the
// lines after the rule under "Search for:", from column 2 to column 77,
// trailing spaces removed, until the first blank line
auto recorded_list(std::string_view const name) -> std::vector<std::string> {
    auto const text = read_file(legacy_dir / "view_edit" / (std::string{name} + ".screens.txt"));
    auto stream = std::istringstream{text.substr(text.rfind("==== "))};
    auto lines = std::vector<std::string>{};
    auto inside = false;
    for (auto line = std::string{}; std::getline(stream, line);) {
        if (line.find("──────────") != std::string::npos and line.find("Search for") == std::string::npos and
            line.starts_with("│")) {
            inside = true;
            continue;
        }
        if (not inside) {
            continue;
        }
        constexpr auto first_byte = std::size_t{4};  // "│" (3 bytes in UTF-8), then a space
        constexpr auto content_columns = std::size_t{76};
        auto content = line.size() > first_byte ? line.substr(first_byte, content_columns) : std::string{};
        content.erase(content.find_last_not_of(' ') + 1);
        if (content.empty()) {
            break;
        }
        lines.push_back(content);
    }
    return lines;
}

auto shown(std::vector<Nucleotide> const& molecule, std::string_view const letters) -> std::vector<std::string> {
    auto const sites = legacy_find_sequence(molecule, letters, Topology::linear);
    return sites.has_value() ? legacy_site_list_lines(*sites) : std::vector<std::string>{};
}

}  // namespace

// Reference: PRIMER.exe, View/Edit -> Sequence -> Search on
// pairs/TEST600.GB (tests/legacy/view_edit/find_*.screens.txt, recorded
// 2026-10-21): both strands (TTAAAT: 8 and 589C), IUPAC letters
// (CATNTA: 5 sites), six sites a line (AAAT: 10 sites)
TEST_CASE("Find Sequence: the recorded site lists", "[search][view][legacy]") {
    auto const molecule = test600();
    for (auto const& [name, letters] : {std::pair{"find_catcta", "CATCTA"}, std::pair{"find_ttaaat", "TTAAAT"},
                                        std::pair{"find_catnta", "CATNTA"}, std::pair{"find_aaat", "AAAT"}}) {
        INFO(name);
        auto const recorded = recorded_list(name);
        REQUIRE_FALSE(recorded.empty());
        CHECK(shown(molecule, letters) == recorded);
    }
}

// Same reference: "Number of Sites =  10" for AAAT
TEST_CASE("Find Sequence: the number of sites", "[search][view][legacy]") {
    auto const sites = legacy_find_sequence(test600(), "aaat", Topology::linear);
    REQUIRE(sites.has_value());
    CHECK(sites.value_or(std::vector<primerpicker::enzymes::Site>{}).size() == 10);
    CHECK_FALSE(legacy_find_sequence(test600(), "AXT", Topology::linear).has_value());
}
