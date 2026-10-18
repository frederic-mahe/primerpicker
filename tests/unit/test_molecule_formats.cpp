#include "core/io/genbank.hpp"
#include "core/io/molecule_file.hpp"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

using primerpicker::io::Molecule;
using primerpicker::io::MoleculeError;
using primerpicker::io::MoleculeResult;
using primerpicker::io::parse_fasta;
using primerpicker::io::parse_molecule_file;
using primerpicker::io::parse_plain_sequence;

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

// throws (and fails the test) when there is no molecule
auto molecule_of(MoleculeResult const& result) -> Molecule {
    if (not result.molecule.has_value()) {
        throw std::runtime_error{"no molecule"};
    }
    return *result.molecule;
}

// first line of TEST600 as PRIMER.exe shows it (Modify Primer Position)
constexpr auto test600_start = std::string_view{"GGTCATTTTAAATGGACAACCATCTATTGAGAATACTCATGTCAACCTGGCTCAAAGGCT"};

}  // namespace

// Reference values: what PRIMER.exe shows after loading each file
// (tests/legacy/README.md, "Molecule file formats")

TEST_CASE("a FASTA molecule takes the first word of its header", "[io][legacy]") {
    auto const result = parse_fasta(legacy_file("molecules/TFASTA.FA"));
    CHECK(result.error == MoleculeError::none);
    auto const molecule = molecule_of(result);
    CHECK(molecule.name == "TEST600");  // PRIMER.exe
    CHECK(molecule.sequence.size() == 600);
    CHECK(molecule.sequence.starts_with(test600_start));  // PRIMER.exe
}

TEST_CASE("long FASTA names are kept, lowercase sequences read", "[io][legacy]") {
    auto const molecule = molecule_of(parse_fasta(legacy_file("molecules/TFA2.FA")));
    CHECK(molecule.name == "VERYLONGMOLECULENAME_2026");  // PRIMER.exe
    CHECK(molecule.sequence.starts_with(test600_start));   // PRIMER.exe
}

TEST_CASE("a FASTA header with spaces gives its first word", "[io][legacy]") {
    CHECK(molecule_of(parse_fasta(legacy_file("molecules/TFA3.FA"))).name == "name");  // PRIMER.exe
}

TEST_CASE("a plain sequence takes the file name", "[io][legacy]") {
    auto const text = legacy_file("molecules/TPLAIN.SEQ");
    auto const molecule = molecule_of(parse_plain_sequence({.text = text, .name = "TPLAIN.SEQ"}));
    CHECK(molecule.name == "TPLAIN.SEQ");  // PRIMER.exe
    CHECK(molecule.sequence.size() == 600);
    CHECK(molecule.sequence.starts_with(test600_start));  // PRIMER.exe
}

TEST_CASE("a FASTA text without header is not FASTA", "[io]") {
    auto const result = parse_fasta("ACGT\n");
    CHECK_FALSE(result.molecule.has_value());
    CHECK(result.error == MoleculeError::missing_header);
}

TEST_CASE("an empty sequence is an error", "[io]") {
    CHECK(parse_fasta(">empty\n").error == MoleculeError::empty_sequence);
    CHECK(parse_plain_sequence({.text = "\r\n", .name = "EMPTY.SEQ"}).error == MoleculeError::empty_sequence);
}

TEST_CASE("the file format is recognised from the content", "[io]") {
    // GenBank, FASTA and plain sequence load the same molecule
    auto const genbank_text = legacy_file("pairs/TEST600.GB");
    auto const fasta_text = legacy_file("molecules/TFASTA.FA");
    auto const plain_text = legacy_file("molecules/TPLAIN.SEQ");
    auto const genbank = molecule_of(parse_molecule_file({.text = genbank_text, .name = "TEST600.GB"}));
    auto const fasta = molecule_of(parse_molecule_file({.text = fasta_text, .name = "TFASTA.FA"}));
    auto const plain = molecule_of(parse_molecule_file({.text = plain_text, .name = "TPLAIN.SEQ"}));
    CHECK(genbank.name == "TEST600");
    CHECK(fasta.name == "TEST600");
    CHECK(plain.name == "TPLAIN.SEQ");
    CHECK(fasta.sequence == genbank.sequence);
    CHECK(plain.sequence == genbank.sequence);
}
