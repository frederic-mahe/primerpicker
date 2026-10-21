#include "core/io/molecule_file.hpp"
#include "core/search/legacy_picking.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using primerpicker::search::legacy_change_length;
using primerpicker::search::PickingMessage;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::PrimerLocation;
using primerpicker::sequence::Strand;

auto test600() -> std::vector<Nucleotide> {
    auto const path = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "pairs" / "TEST600.GB";
    auto file = std::ifstream{path, std::ios::binary};
    auto text = std::ostringstream{};
    text << file.rdbuf();
    auto const parsed = primerpicker::io::parse_molecule_file({.text = std::move(text).str(), .name = "TEST600.GB"});
    auto sequence = parsed.molecule.has_value() ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                                                : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot parse TEST600.GB"};
    }
    return std::move(*sequence);
}

auto letters(std::vector<Nucleotide> const& molecule, PrimerLocation const& location) -> std::string {
    auto const primer = primerpicker::sequence::primer_at(molecule, location);
    auto text = std::string{};
    if (primer.has_value()) {
        std::ranges::transform(*primer, std::back_inserter(text), primerpicker::sequence::to_letter);
    }
    return text;
}

constexpr auto top_51 = PrimerLocation{.position = 51, .length = 20, .strand = Strand::top};

}  // namespace

// Reference: PRIMER.exe, Modify -> Length on the primer at 51 of
// TEST600 (tests/legacy/modify_length/, recorded 2026-10-21): 25 gives
// CTCAAAGGCTCAGACATCAGGGCTC, 1 gives C (the 5' end stays)
TEST_CASE("Modify -> Length keeps the primer's 5' end", "[search][legacy]") {
    auto const molecule = test600();
    auto const longer = legacy_change_length(molecule, top_51, 25);
    CHECK(longer.message == PickingMessage::none);
    CHECK(letters(molecule, longer.location) == "CTCAAAGGCTCAGACATCAGGGCTC");
    CHECK(letters(molecule, legacy_change_length(molecule, top_51, 1).location) == "C");
}

// Same reference: 400(C), 25 bases: TGTATCCAGGTAATGCAAGGATACC
TEST_CASE("Modify -> Length on the complement strand", "[search][legacy]") {
    auto const molecule = test600();
    constexpr auto complement_400 = PrimerLocation{.position = 400, .length = 20, .strand = Strand::complement};
    auto const changed = legacy_change_length(molecule, complement_400, 25);
    CHECK(changed.location.position == 400);
    CHECK(letters(molecule, changed.location) == "TGTATCCAGGTAATGCAAGGATACC");
}

// Same reference: 0 gives "Primer must have length!"; 30 bases at 575
// "Primer Too Long - Extends outside molecule"; the primer is kept
TEST_CASE("Modify -> Length refuses no base and a primer off the molecule", "[search][legacy]") {
    auto const molecule = test600();
    auto const none = legacy_change_length(molecule, top_51, 0);
    CHECK(none.message == PickingMessage::primer_must_have_length);
    CHECK(none.location.length == 20);
    constexpr auto top_575 = PrimerLocation{.position = 575, .length = 20, .strand = Strand::top};
    auto const too_long = legacy_change_length(molecule, top_575, 30);
    CHECK(too_long.message == PickingMessage::primer_too_long);
    CHECK(too_long.location.length == 20);
}
