#include "core/analysis/legacy_verdict.hpp"
#include "core/io/molecule_file.hpp"
#include "core/search/legacy_picking.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using primerpicker::analysis::LegacyCriteria;
using primerpicker::search::legacy_find_next;
using primerpicker::search::legacy_picking_move;
using primerpicker::search::legacy_picking_state;
using primerpicker::search::PickingKey;
using primerpicker::search::PickingMessage;
using primerpicker::search::PickingState;
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

auto molecule(std::filesystem::path const& relative_path) -> std::vector<Nucleotide> {
    auto const text = read_file(legacy_dir / relative_path);
    auto const parsed =
        primerpicker::io::parse_molecule_file({.text = text, .name = relative_path.filename().string()});
    auto sequence = parsed.molecule.has_value() ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                                                : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot parse " + relative_path.string()};
    }
    return std::move(*sequence);
}

auto test600() -> std::vector<Nucleotide> {
    return molecule("pairs/TEST600.GB");
}

constexpr auto top(std::size_t const position, std::size_t const length = 20) -> PrimerLocation {
    return {.position = position, .length = length, .strand = Strand::top};
}

constexpr auto complement(std::size_t const position, std::size_t const length = 20) -> PrimerLocation {
    return {.position = position, .length = length, .strand = Strand::complement};
}

auto same(PrimerLocation const& one, PrimerLocation const& other) -> bool {
    return one.position == other.position and one.length == other.length and one.strand == other.strand;
}

// A recorded state of the picking screen
struct RecordedState {
    std::string run;
    std::string molecule_file;
    PrimerLocation location;
    bool meets;
    std::string value_line;
    std::string window_colour;
};

auto split(std::string const& line) -> std::vector<std::string> {
    auto fields = std::vector<std::string>{};
    auto stream = std::istringstream{line};
    for (auto field = std::string{}; std::getline(stream, field, '\t');) {
        fields.push_back(field);
    }
    return fields;
}

// tests/legacy/picking/picking.tsv: run, molecule file (from
// tests/legacy), frame, strand, length, position, meets, value line,
// window colour
auto recorded_states() -> std::vector<RecordedState> {
    auto file = std::ifstream{legacy_dir / "picking" / "picking.tsv"};
    auto states = std::vector<RecordedState>{};
    constexpr auto field_count = std::size_t{9};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto const fields = split(line);
        if (fields.size() != field_count) {
            throw std::runtime_error{"cannot read line: " + line};
        }
        auto const strand = fields[3] == "complement" ? Strand::complement : Strand::top;
        states.push_back({.run = fields[0],
                          .molecule_file = fields[1],
                          .location = {.position = std::stoul(fields[5]),
                                       .length = std::stoul(fields[4]),
                                       .strand = strand},
                          .meets = fields[6] == "Yes",
                          .value_line = fields[7],
                          .window_colour = fields[8]});
    }
    return states;
}

void check_recorded_state(RecordedState const& recorded) {
    INFO(recorded.run << ": " << recorded.value_line);
    auto const sequence = molecule(recorded.molecule_file);
    auto const shown = legacy_picking_state(sequence, recorded.location, LegacyCriteria{});
    REQUIRE(shown.has_value());
    CHECK(shown.value_or(PickingState{}).value_line == recorded.value_line);
    CHECK(shown.value_or(PickingState{}).meets_criteria == recorded.meets);
    CHECK((recorded.window_colour == "green") == recorded.meets);
}

// The position found by Find next, 0 for "No Primers found"
auto found_position(std::optional<PrimerLocation> const& found) -> std::size_t {
    return found.has_value() ? found->position : std::size_t{0};
}

auto found_length(std::optional<PrimerLocation> const& found) -> std::size_t {
    return found.has_value() ? found->length : std::size_t{0};
}

auto found_strand(std::optional<PrimerLocation> const& found) -> std::optional<Strand> {
    return found.has_value() ? std::optional{found->strand} : std::nullopt;
}

}  // namespace

// Reference: PRIMER.exe, From Molecule (built-in criteria), every
// settled state of the picking screen recorded on 2026-10-21
// (tests/legacy/picking/picking.tsv, tools/picking_table.py): the value
// line, MEETS CRITERIA SELECTED, and the window's background (green
// when the primer meets the criteria, light grey otherwise)
TEST_CASE("the picking screen's values, verdict and window colour, as recorded", "[search][picking][legacy]") {
    auto const states = recorded_states();
    REQUIRE(states.size() >= 40);
    for (auto const& recorded : states) {
        check_recorded_state(recorded);
    }
}

// Reference: PRIMER.exe, TEST600, primer at 51 (20 bases):
// tests/legacy/picking/explore_51.screens.txt (right right right left:
// 52, 53, 54, 53; Down from 53: 113, then Up: 53)
TEST_CASE("the arrows move the primer by one base, Up and Down by a line of 60", "[search][picking][legacy]") {
    constexpr auto size = std::size_t{600};
    CHECK(same(legacy_picking_move(top(51), size, PickingKey::right).location, top(52)));
    CHECK(same(legacy_picking_move(top(54), size, PickingKey::left).location, top(53)));
    CHECK(same(legacy_picking_move(top(53), size, PickingKey::down).location, top(113)));
    CHECK(same(legacy_picking_move(top(113), size, PickingKey::up).location, top(53)));
    CHECK(legacy_picking_move(top(51), size, PickingKey::right).message == PickingMessage::none);
}

// Reference: PRIMER.exe, TEST600: from 2, Left gives 1 and a second
// Left keeps it; Up from 51 does nothing (start_2, pages_51); from 579,
// Right gives 580, 581, then nothing; Down from 577 does nothing
// (end_579, end_575); no message
TEST_CASE("a move that would take the primer off the molecule is ignored", "[search][picking][legacy]") {
    constexpr auto size = std::size_t{600};
    CHECK(same(legacy_picking_move(top(2), size, PickingKey::left).location, top(1)));
    CHECK(same(legacy_picking_move(top(1), size, PickingKey::left).location, top(1)));
    CHECK(same(legacy_picking_move(top(51), size, PickingKey::up).location, top(51)));
    CHECK(same(legacy_picking_move(top(580), size, PickingKey::right).location, top(581)));
    CHECK(same(legacy_picking_move(top(581), size, PickingKey::right).location, top(581)));
    CHECK(same(legacy_picking_move(top(577), size, PickingKey::down).location, top(577)));
    CHECK(legacy_picking_move(top(581), size, PickingKey::right).message == PickingMessage::none);
    CHECK(legacy_picking_move(top(1), size, PickingKey::left).message == PickingMessage::none);
}

// Reference: PRIMER.exe, TEST600, complement strand at 400(C) (comp_400):
// Right gives 401, Left 400, Down 460: positions grow to the right on
// both strands
TEST_CASE("on the complement strand, Right and Down raise the position too", "[search][picking][legacy]") {
    constexpr auto size = std::size_t{600};
    CHECK(same(legacy_picking_move(complement(400), size, PickingKey::right).location, complement(401)));
    CHECK(same(legacy_picking_move(complement(401), size, PickingKey::left).location, complement(400)));
    CHECK(same(legacy_picking_move(complement(400), size, PickingKey::down).location, complement(460)));
    CHECK(same(legacy_picking_move(complement(20), size, PickingKey::left).location, complement(20)));
    CHECK(same(legacy_picking_move(complement(600), size, PickingKey::right).location, complement(600)));
}

// Reference: PRIMER.exe, TEST600 at 53 then 36: + and - change the
// length by one base (explore_51, len_plus); + at 581 (20 bases) shows
// "Primer Too Long - Extends outside molecule" and keeps the primer
// (end_579)
TEST_CASE("+ and - change the primer's length", "[search][picking][legacy]") {
    constexpr auto size = std::size_t{600};
    CHECK(same(legacy_picking_move(top(53), size, PickingKey::longer).location, top(53, 21)));
    CHECK(same(legacy_picking_move(top(53, 21), size, PickingKey::shorter).location, top(53, 20)));
    auto const too_long = legacy_picking_move(top(581), size, PickingKey::longer);
    CHECK(same(too_long.location, top(581)));
    CHECK(too_long.message == PickingMessage::primer_too_long);
    auto const too_long_complement = legacy_picking_move(complement(20), size, PickingKey::longer);
    CHECK(same(too_long_complement.location, complement(20)));
    CHECK(too_long_complement.message == PickingMessage::primer_too_long);
}

// Reference: PRIMER.exe, TEST600, built-in criteria (find_51, find_end,
// len_plus, none_2): Ctrl+Right from 51 gives 71, then 72; Ctrl+Left
// from 72 gives 71, 36, 35; from 530, 559 then "No Primers found"; with
// 22 bases, Ctrl+Right from 36 gives 71; from 2, Ctrl+Left finds none
TEST_CASE("Ctrl+arrows find the next primer that meets the criteria", "[search][picking][legacy]") {
    auto const sequence = test600();
    auto const criteria = LegacyCriteria{};
    auto const next = [&](PrimerLocation const& from, SearchDirection const direction) {
        return found_position(legacy_find_next(sequence, from, criteria, direction));
    };
    CHECK(next(top(51), SearchDirection::forward) == 71);
    CHECK(next(top(71), SearchDirection::forward) == 72);
    CHECK(next(top(72), SearchDirection::backward) == 71);
    CHECK(next(top(71), SearchDirection::backward) == 36);
    CHECK(next(top(36), SearchDirection::backward) == 35);
}

// Same references
TEST_CASE("Ctrl+arrows: no primer found, and the length is kept", "[search][picking][legacy]") {
    auto const sequence = test600();
    auto const criteria = LegacyCriteria{};
    auto const next = [&](PrimerLocation const& from, SearchDirection const direction) {
        return found_position(legacy_find_next(sequence, from, criteria, direction));
    };
    CHECK(next(top(530), SearchDirection::forward) == 559);
    CHECK(next(top(559), SearchDirection::forward) == 0);
    CHECK(next(top(36, 22), SearchDirection::forward) == 71);
    CHECK(found_length(legacy_find_next(sequence, top(36, 22), criteria, SearchDirection::forward)) == 22);
    CHECK(next(top(2), SearchDirection::backward) == 0);
}

// Reference: PRIMER.exe, TEST600, complement strand (comp_400): from
// 460(C), Ctrl+Right gives 516(C); from 516(C), Ctrl+Left gives 457(C)
TEST_CASE("Ctrl+arrows on the complement strand", "[search][picking][legacy]") {
    auto const sequence = test600();
    auto const criteria = LegacyCriteria{};
    auto const forward = legacy_find_next(sequence, complement(460), criteria, SearchDirection::forward);
    auto const backward = legacy_find_next(sequence, complement(516), criteria, SearchDirection::backward);
    CHECK(found_position(forward) == 516);
    CHECK(found_strand(forward) == Strand::complement);
    CHECK(found_position(backward) == 457);
    CHECK(found_strand(backward) == Strand::complement);
}

// Reference: PRIMER.exe, TEST600 and picking/TWIDTH.GB: - on a primer of
// one base shows "Primer must have length!" (w_short_41: 5, 4, 3, 2, 1
// bases, then the message)
TEST_CASE("the - key stops at one base", "[search][picking][legacy]") {
    constexpr auto size = std::size_t{300};
    CHECK(same(legacy_picking_move(top(41, 2), size, PickingKey::shorter).location, top(41, 1)));
    auto const refused = legacy_picking_move(top(41, 1), size, PickingKey::shorter);
    CHECK(same(refused.location, top(41, 1)));
    CHECK(refused.message == PickingMessage::primer_must_have_length);
}

// Reference: PRIMER.exe, search/TRUN.GB (trun_450, trun_find_440,
// recorded 2026-10-21): the primer at 450 meets every criterion of the
// Primer Analysis screen but ends with AAA: "Runs: !", MEETS CRITERIA
// SELECTED: No, a light grey window; Ctrl+Right from 449 skips it (and
// 451, AAA too) and gives 471; from 40, 42 then 46 (trun_find_40)
TEST_CASE("the run at the 3' end counts on the picking screen", "[search][picking][legacy]") {
    auto const sequence = molecule("search/TRUN.GB");
    auto const criteria = LegacyCriteria{};
    auto const state = legacy_picking_state(sequence, top(450), criteria).value_or(PickingState{});
    CHECK_FALSE(state.meets_criteria);
    CHECK(state.value_line.ends_with("Runs: !"));
    CHECK(found_position(legacy_find_next(sequence, top(449), criteria, SearchDirection::forward)) == 471);
    CHECK(found_position(legacy_find_next(sequence, top(40), criteria, SearchDirection::forward)) == 42);
}
