#include "core/analysis/legacy_verdict.hpp"
#include "core/io/molecule_file.hpp"
#include "core/search/legacy_picking.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

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
#include <utility>
#include <vector>

namespace {

using primerpicker::analysis::LegacyCriteria;
using primerpicker::search::legacy_find_next;
using primerpicker::search::legacy_picking_state;
using primerpicker::search::SearchDirection;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::PrimerLocation;
using primerpicker::sequence::Strand;

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

auto nucleotides(std::string const& letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + letters};
    }
    return std::move(*converted);
}

auto test600() -> std::vector<Nucleotide> {
    auto const text = read_file(legacy_dir / "pairs" / "TEST600.GB");
    auto const parsed = primerpicker::io::parse_molecule_file({.text = text, .name = "TEST600.GB"});
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"cannot parse TEST600.GB"};
    }
    return nucleotides(parsed.molecule->sequence);
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

// A recorded state of the picking screen of a primer with a 5' extension
struct RecordedState {
    std::size_t length;  // extension included, as "Primer length"
    std::vector<Nucleotide> extension;
    std::size_t position;
    bool meets;
    std::string value_line;
};

// extension/picking_extension.tsv: run, molecule, frame, strand, length,
// extension, position, meets, value line, window colour
auto recorded_states() -> std::vector<RecordedState> {
    constexpr auto column_count = 10U;
    auto states = std::vector<RecordedState>{};
    auto file = std::ifstream{legacy_dir / "extension" / "picking_extension.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        if (fields.size() != column_count or fields[1] != "pairs/TEST600.GB" or fields[3] != "top") {
            throw std::runtime_error{"unexpected line: " + line};
        }
        states.push_back({.length = number(fields[4]),
                          .extension = nucleotides(fields[5]),
                          .position = number(fields[6]),
                          .meets = fields[7] == "Yes",
                          .value_line = fields[8]});
    }
    return states;
}

}  // namespace

// Reference values: PRIMER.exe, the picking screen of a primer made with
// From Molecule at 51 of TEST600 with the 5' extension GAATTC, 20 bases
// in all (tests/legacy/extension/picking_extension.tsv, recorded
// 2026-10-22): the values are those of the whole primer, extension
// included (the Tm without it), the window covers the molecule's bases

TEST_CASE("picking screen, 5' extension: every recorded state", "[search][legacy]") {
    auto const molecule = test600();
    auto const states = recorded_states();
    REQUIRE(states.size() == 5);
    for (auto const& state : states) {
        CAPTURE(state.position, state.length);
        auto const location = PrimerLocation{
            .position = state.position, .length = state.length - state.extension.size(), .strand = Strand::top};
        auto const shown = legacy_picking_state(molecule, location, LegacyCriteria{}, state.extension)
                               .value_or(primerpicker::search::PickingState{.value_line = "none", .meets_criteria = false});
        CHECK(shown.value_line == state.value_line);
        CHECK(shown.meets_criteria == state.meets);
    }
}

TEST_CASE("picking screen, 5' extension: Ctrl+Right finds a whole primer that meets the criteria",
          "[search][legacy]") {
    // PRIMER.exe (extension/ext_picking3): from 52, 14 bases after
    // GAATTC, Ctrl+Right moves to 75
    auto const molecule = test600();
    auto const extension = nucleotides("GAATTC");
    auto const found = legacy_find_next(molecule, {.position = 52, .length = 14, .strand = Strand::top},
                                        LegacyCriteria{}, SearchDirection::forward, extension);
    CHECK(found.value_or(PrimerLocation{.position = 0, .length = 0, .strand = Strand::top}).position == 75);
    CHECK(found.value_or(PrimerLocation{.position = 0, .length = 0, .strand = Strand::top}).length == 14);
}
