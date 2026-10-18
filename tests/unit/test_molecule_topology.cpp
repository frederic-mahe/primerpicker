#include "core/io/genbank.hpp"
#include "core/io/molecule_file.hpp"
#include "core/sequence/topology.hpp"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

using primerpicker::io::parse_genbank;
using primerpicker::io::parse_molecule_file;
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

auto genbank_topology(std::string_view const text) -> Topology {
    auto const parsed = parse_genbank(text);
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"not a GenBank text"};
    }
    return parsed.molecule->topology;
}

}  // namespace

// Reference: PRIMER.exe searches TENZC.GB ("circular" on its LOCUS
// line) across the origin, and TENZ.GB ("linear") not
// (tests/legacy/README.md, "True and almost sites")

TEST_CASE("GenBank: the LOCUS line gives the topology", "[io][genbank]") {
    auto const legacy = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "enzymes";
    CHECK(genbank_topology(read_file(legacy / "TENZC.GB")) == Topology::circular);
    CHECK(genbank_topology(read_file(legacy / "TENZ.GB")) == Topology::linear);
}

TEST_CASE("GenBank: without a topology on the LOCUS line, a molecule is linear", "[io][genbank]") {
    CHECK(genbank_topology("LOCUS       ABC\nORIGIN\n        1 acgt\n//\n") == Topology::linear);
}

TEST_CASE("GenBank: \"circular\" is read in either case", "[io][genbank]") {
    // not recorded: PRIMER.exe was only given lower case
    CHECK(genbank_topology("LOCUS       ABC  4 bp  DNA  CIRCULAR\nORIGIN\n        1 acgt\n//\n") ==
          Topology::circular);
}

TEST_CASE("GenBank: \"circular\" counts on the LOCUS line only", "[io][genbank]") {
    CHECK(genbank_topology("LOCUS       ABC  4 bp  DNA  linear\nDEFINITION  a circular plasmid\nORIGIN\n"
                           "        1 acgt\n//\n") == Topology::linear);
}

TEST_CASE("FASTA and plain sequences are linear", "[io]") {
    auto const fasta = parse_molecule_file({.text = ">ABC\nACGT\n", .name = "ABC.FA"});
    REQUIRE(fasta.molecule.has_value());
    CHECK(fasta.molecule.value_or(primerpicker::io::Molecule{}).topology == Topology::linear);
    auto const plain = parse_molecule_file({.text = "ACGT\n", .name = "ABC.SEQ"});
    REQUIRE(plain.molecule.has_value());
    CHECK(plain.molecule.value_or(primerpicker::io::Molecule{}).topology == Topology::linear);
}
