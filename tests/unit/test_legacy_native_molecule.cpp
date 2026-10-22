#include "core/io/genbank.hpp"
#include "core/io/legacy_native_molecule.hpp"
#include "core/io/molecule_file.hpp"
#include "core/sequence/topology.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

using primerpicker::io::ComplementFlag;
using primerpicker::io::GeneDirection;
using primerpicker::io::MarkerDirection;
using primerpicker::io::Molecule;
using primerpicker::io::MoleculeError;
using primerpicker::io::NativeEnds;
using primerpicker::io::NativeError;
using primerpicker::io::NativeGene;
using primerpicker::io::NativeMolecule;
using primerpicker::io::SequenceKind;
using primerpicker::io::native_from_molecule;
using primerpicker::io::parse_molecule_file;
using primerpicker::io::parse_native_molecule;
using primerpicker::io::shows_complement;
using primerpicker::io::write_native_molecule;
using primerpicker::sequence::Topology;

auto read_file(std::filesystem::path const& path) -> std::string {
    auto text = std::string(std::filesystem::file_size(path), '\0');
    auto file = std::ifstream{path, std::ios::binary};
    file.read(text.data(), std::ssize(text));
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    return text;
}

auto legacy_file(std::string_view const relative) -> std::string {
    return read_file(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / relative);
}

auto native(std::string_view const relative) -> NativeMolecule {
    auto const result = parse_native_molecule(legacy_file(relative));
    if (not result.molecule.has_value()) {
        throw std::runtime_error{"not read: " + std::string{relative}};
    }
    return *result.molecule;
}

// throws (and fails the test) when the molecule is circular
auto ends_of(NativeMolecule const& molecule) -> NativeEnds {
    if (not molecule.linear_ends.has_value()) {
        throw std::runtime_error{"no linear ends: " + molecule.name};
    }
    return *molecule.linear_ends;
}

// the text written back from a file read; throws when it is not read
auto written_back(std::string_view const text) -> std::string {
    auto const result = parse_native_molecule(text);
    if (not result.molecule.has_value()) {
        throw std::runtime_error{"not read"};
    }
    return write_native_molecule(*result.molecule);
}

auto loaded(std::string_view const relative) -> Molecule {
    auto const path = std::filesystem::path{relative};
    auto const text = legacy_file(relative);
    auto const name = path.filename().string();
    auto const result = parse_molecule_file({.text = text, .name = name});
    if (not result.molecule.has_value()) {
        throw std::runtime_error{"not loaded: " + std::string{relative}};
    }
    return *result.molecule;
}

// the files PRIMER.exe wrote (Save/Remove -> S, Create New)
constexpr auto written_files = std::to_array<std::string_view>(
    {"molecule_save/written/save_new.ZZ.GB", "molecule_save/written/save_fa.TFASTA.FA",
     "molecule_save/written/save_circ.ZC.GB", "molecule_save/written/create_nothing.ZM",
     "molecule_save/written/create_feat4.ZN", "molecule_save/written/create_prot.ZP",
     "molecule_save/written/feat_marker4.ZF", "molecule_save/written/feat_region2.ZF",
     "molecule_save/written/feat_gback.ZF"});

}  // namespace

// Reference values: PRIMER.exe, Molecule-Ops -> Save/Remove (S) and
// Create New (tests/legacy/molecule_save/, recorded 2026-10-22;
// tests/legacy/PROBES_20261022_molecule_save.md)

TEST_CASE("native molecules: every file PRIMER.exe wrote is read and written back byte for byte",
          "[io][legacy]") {
    for (auto const relative : written_files) {
        CAPTURE(relative);
        auto const text = legacy_file(relative);
        CHECK(parse_native_molecule(text).error == NativeError::none);
        CHECK(written_back(text) == text);
    }
}

TEST_CASE("native molecules: Save writes a GenBank or FASTA molecule as PRIMER.exe does", "[io][legacy]") {
    // save_new: TEST600.GB saved as ZZ.GB; save_fa: TFASTA.FA overwritten;
    // save_circ: TCIRC.GB (circular) saved as ZC.GB
    CHECK(write_native_molecule(native_from_molecule(loaded("pairs/TEST600.GB"))) ==
          legacy_file("molecule_save/written/save_new.ZZ.GB"));
    CHECK(write_native_molecule(native_from_molecule(loaded("molecules/TFASTA.FA"))) ==
          legacy_file("molecule_save/written/save_fa.TFASTA.FA"));
    CHECK(write_native_molecule(native_from_molecule(loaded("circular/TCIRC.GB"))) ==
          legacy_file("molecule_save/written/save_circ.ZC.GB"));
}

TEST_CASE("native molecules: the fields of a created molecule", "[io][legacy]") {
    // create_feat4: name zn, size 12, left end e, right end b, description
    // de, gene go (description do) from 5 to 9, no sequence
    auto const created = native("molecule_save/written/create_feat4.ZN");
    CHECK(created.name == "zn");
    CHECK(created.size == 12);
    REQUIRE(created.genes.size() == 1);
    CHECK(created.genes.front().name == "go");
    CHECK(created.genes.front().description == "do");
    CHECK(created.genes.front().start == 5);
    CHECK(created.genes.front().end == 9);
    CHECK(created.genes.front().direction == GeneDirection::clockwise);
    CHECK(created.genes.front().flag == ComplementFlag::no);
    CHECK(created.markers.empty());
    CHECK(ends_of(created).left == "e");
    CHECK(ends_of(created).right == "b");
    CHECK(created.description == "de");
    CHECK_FALSE(created.sequence.has_value());
}

TEST_CASE("native molecules: markers and circular molecules", "[io][legacy]") {
    // feat_marker4: marker go at 5, Complement Y
    auto const marked = native("molecule_save/written/feat_marker4.ZF");
    REQUIRE(marked.markers.size() == 1);
    CHECK(marked.markers.front().position == 5);
    CHECK(marked.markers.front().direction == MarkerDirection::left);
    // save_circ: circular, no linear block
    auto const circular = native("molecule_save/written/save_circ.ZC.GB");
    CHECK_FALSE(circular.linear_ends.has_value());
    CHECK(circular.sequence.has_value());
}

TEST_CASE("native molecules: files PRIMER.exe reads", "[io][legacy]") {
    // hand7 (inputs/ZH7): a gene and two markers, a description and a
    // sequence; PRIMER.exe lists GENE 2 6 (C), MARKER 3 », MARKER 8 « (C)
    auto const both = native("molecule_save/inputs/ZH7");
    CHECK(both.genes.size() == 1);
    CHECK(both.markers.size() == 2);
    CHECK(both.description == "both");
    CHECK(both.sequence == "ACGTACGTGGCC");
    // hand3 (inputs/ZH3): circular; a gene line without NUL is shown as
    // its name and its description (gcdc, gcdc); a flag R loads
    auto const circular = native("molecule_save/inputs/ZH3");
    CHECK_FALSE(circular.linear_ends.has_value());
    REQUIRE(circular.genes.size() == 1);
    CHECK(circular.genes.front().name == "gcdc");
    CHECK(circular.genes.front().description == "gcdc");
    // hand4 (inputs/ZH4): sequenceAA, "4 AA PROTEIN"
    CHECK(native("molecule_save/inputs/ZH4").kind == SequenceKind::protein);
}

TEST_CASE("native molecules: files PRIMER.exe rejects (Illegal File Format)", "[io][legacy]") {
    // hand1, hand9, handA: a third line other than 0; hand2: marker
    // records that are not markers
    for (auto const* const relative :
         {"molecule_save/inputs/ZH1", "molecule_save/inputs/ZH2", "molecule_save/inputs/ZH9",
          "molecule_save/inputs/ZHA"}) {
        CAPTURE(relative);
        auto const result = parse_native_molecule(legacy_file(relative));
        CHECK_FALSE(result.molecule.has_value());
        CHECK(result.error == NativeError::illegal_format);
    }
}

TEST_CASE("native molecules: when the features table shows (C) for a gene", "[io][legacy]") {
    // PRIMER.exe, View/Edit -> Features (hand6, hand7, hand8, feat_gback,
    // reload_zn): (C) when the flag is Y and start < end, or N and
    // start > end; cw or ccw does not matter. hand8:
    auto const gene = [](std::size_t const start, std::size_t const end, GeneDirection const direction,
                         ComplementFlag const flag) {
        return NativeGene{.name = "g", .description = "d", .start = start, .end = end, .direction = direction, .flag = flag};
    };
    CHECK(shows_complement(gene(2, 6, GeneDirection::clockwise, ComplementFlag::yes)));          // hand8 ga
    CHECK_FALSE(shows_complement(gene(9, 5, GeneDirection::clockwise, ComplementFlag::yes)));    // hand8 gb
    CHECK(shows_complement(gene(9, 5, GeneDirection::counterclockwise, ComplementFlag::no)));    // hand8 gc
    CHECK_FALSE(shows_complement(gene(9, 5, GeneDirection::counterclockwise, ComplementFlag::yes)));  // hand8 gd
}

TEST_CASE("native molecules: genes shown with or without (C)", "[io][legacy]") {
    // PRIMER.exe, View/Edit -> Features: hand7, hand6, reload_zn (as above)
    auto const gene = [](std::size_t const start, std::size_t const end, GeneDirection const direction,
                         ComplementFlag const flag) {
        return NativeGene{.name = "g", .description = "d", .start = start, .end = end, .direction = direction, .flag = flag};
    };
    CHECK(shows_complement(gene(2, 6, GeneDirection::counterclockwise, ComplementFlag::yes)));   // hand7 ga
    CHECK_FALSE(shows_complement(gene(2, 6, GeneDirection::counterclockwise, ComplementFlag::no)));  // hand6 ga
    CHECK_FALSE(shows_complement(gene(7, 11, GeneDirection::clockwise, ComplementFlag::no)));    // hand6 gb
    CHECK_FALSE(shows_complement(gene(5, 9, GeneDirection::clockwise, ComplementFlag::no)));     // reload_zn
}

TEST_CASE("native molecules: loaded as molecules", "[io][legacy]") {
    // reload_header, reload_primer: ZZ.GB (written by PRIMER.exe) loads as
    // TEST600, 600 bases, linear; hand3: circular ("CIRCULAR")
    auto const reloaded = loaded("molecule_save/written/save_new.ZZ.GB");
    CHECK(reloaded.name == "TEST600");
    CHECK(reloaded.sequence == loaded("pairs/TEST600.GB").sequence);
    CHECK(reloaded.topology == Topology::linear);
    CHECK(loaded("molecule_save/inputs/ZH3").topology == Topology::circular);
    // a protein or a molecule without sequence cannot make primers
    auto const protein = legacy_file("molecule_save/inputs/ZH4");
    CHECK(parse_molecule_file({.text = protein, .name = "ZH4"}).error == MoleculeError::protein);
    auto const empty = legacy_file("molecule_save/written/create_nothing.ZM");
    CHECK(parse_molecule_file({.text = empty, .name = "ZM"}).error == MoleculeError::empty_sequence);
    auto const illegal = legacy_file("molecule_save/inputs/ZH9");
    CHECK(parse_molecule_file({.text = illegal, .name = "ZH9"}).error == MoleculeError::illegal_format);
}
