#include "core/io/genbank.hpp"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

using primerpicker::io::GenbankError;
using primerpicker::io::GenbankResult;
using primerpicker::io::Molecule;
using primerpicker::io::parse_genbank;

// throws (and fails the test) when parsing failed
auto molecule_of(GenbankResult const& result) -> Molecule {
    if (not result.molecule.has_value()) {
        throw std::runtime_error{"no molecule"};
    }
    return *result.molecule;
}

auto read_file(std::filesystem::path const& path) -> std::string {
    auto file = std::ifstream{path, std::ios::binary};
    return {std::istreambuf_iterator<char>{file}, std::istreambuf_iterator<char>{}};
}

auto legacy_molecule(std::string_view const name) -> std::string {
    return read_file(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "pairs" / name);
}

}  // namespace

// Reference values: what PRIMER.exe shows after loading the file
// (tests/legacy/pairs/, recorded 2026-10-17): the molecule name ("From
// molecule: TEST600"), its size ("600 bps") and its first line in the
// Modify Primer Position screen.

TEST_CASE("a GenBank molecule loads with its name and sequence", "[io]") {
    auto const result = parse_genbank(legacy_molecule("TEST600.GB"));
    CHECK(result.error == GenbankError::none);
    auto const molecule = molecule_of(result);
    CHECK(molecule.name == "TEST600");          // PRIMER.exe
    CHECK(molecule.sequence.size() == 600);     // PRIMER.exe
    CHECK(molecule.sequence.starts_with(        // PRIMER.exe
        "GGTCATTTTAAATGGACAACCATCTATTGAGAATACTCATGTCAACCTGGCTCAAAGGCT"));
}

TEST_CASE("the molecule name comes from LOCUS, not from the file name", "[io]") {
    auto const molecule = molecule_of(parse_genbank(legacy_molecule("TGRAD.GB")));
    CHECK(molecule.name == "GRAD");  // PRIMER.exe: "From molecule: GRAD"
    CHECK(molecule.sequence.size() == 1200);  // PRIMER.exe: map scale to 1200
}

TEST_CASE("Unix line ends are accepted", "[io]") {
    constexpr auto text = std::string_view{"LOCUS       TINY 12 bp DNA linear\n"
                                           "ORIGIN\n"
                                           "        1 acgtacgtac gt\n"
                                           "//\n"};
    auto const molecule = molecule_of(parse_genbank(text));
    CHECK(molecule.name == "TINY");
    CHECK(molecule.sequence == "ACGTACGTACGT");
}

TEST_CASE("a text without LOCUS is not a GenBank molecule", "[io]") {
    auto const result = parse_genbank(">fasta\nACGT\n");
    CHECK_FALSE(result.molecule.has_value());
    CHECK(result.error == GenbankError::missing_locus);
}

TEST_CASE("a GenBank molecule needs an ORIGIN section", "[io]") {
    auto const result = parse_genbank("LOCUS       TINY 12 bp DNA linear\n//\n");
    CHECK_FALSE(result.molecule.has_value());
    CHECK(result.error == GenbankError::missing_origin);
}

TEST_CASE("GenBank parsing works at compile time", "[io][constexpr]") {
    STATIC_CHECK(parse_genbank("LOCUS       X 4 bp\nORIGIN\n        1 acgt\n//\n").molecule->sequence.size() == 4);
}
