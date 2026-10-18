#include "core/sequence/iupac.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <optional>

namespace {

using primerpicker::sequence::from_iupac;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::NucleotideSet;

using enum Nucleotide;

auto code(char const letter) -> NucleotideSet {
    auto const set = from_iupac(letter);
    REQUIRE(set.has_value());
    return set.value_or(NucleotideSet{});
}

}  // namespace

// Reference: Nomenclature Committee of the International Union of
// Biochemistry (NC-IUB) 1985, Nomenclature for incompletely specified
// bases in nucleic acid sequences, Eur J Biochem 150:1-5, Table 1

TEST_CASE("IUPAC letters stand for sets of nucleotides", "[sequence][iupac]") {
    CHECK(code('A') == NucleotideSet{adenine});
    CHECK(code('R') == NucleotideSet{adenine, guanine});
    CHECK(code('Y') == NucleotideSet{cytosine, thymine});
    CHECK(code('S') == NucleotideSet{cytosine, guanine});
    CHECK(code('W') == NucleotideSet{adenine, thymine});
    CHECK(code('K') == NucleotideSet{guanine, thymine});
    CHECK(code('M') == NucleotideSet{adenine, cytosine});
    CHECK(code('B') == NucleotideSet{cytosine, guanine, thymine});
    CHECK(code('D') == NucleotideSet{adenine, guanine, thymine});
    CHECK(code('H') == NucleotideSet{adenine, cytosine, thymine});
    CHECK(code('V') == NucleotideSet{adenine, cytosine, guanine});
    CHECK(code('N') == NucleotideSet{adenine, cytosine, guanine, thymine});
    CHECK(code('n') == NucleotideSet{adenine, cytosine, guanine, thymine});  // either case
}

TEST_CASE("Other letters are not IUPAC nucleotide codes", "[sequence][iupac]") {
    CHECK(not from_iupac('U').has_value());
    CHECK(not from_iupac('X').has_value());
    CHECK(not from_iupac('-').has_value());
}

TEST_CASE("IUPAC sets: membership and complement", "[sequence][iupac]") {
    CHECK(code('R').contains(adenine));
    CHECK(code('R').contains(guanine));
    CHECK(not code('R').contains(cytosine));
    // complement of R (A or G) is Y (T or C); of B (not A) is V (not T)
    CHECK(code('R').complement() == code('Y'));
    CHECK(code('M').complement() == code('K'));
    CHECK(code('B').complement() == code('V'));
    CHECK(code('S').complement() == code('S'));
    CHECK(code('N').complement() == code('N'));
}

TEST_CASE("IUPAC sequences: reverse complement", "[sequence][iupac]") {
    using primerpicker::sequence::reverse_complement_iupac;
    CHECK(reverse_complement_iupac("GGATC") == "GATCC");
    CHECK(reverse_complement_iupac("GTMKAC") == "GTMKAC");  // M (A, C) pairs with K (G, T)
    CHECK(reverse_complement_iupac("ACNNNNGTAYC") == "GRTACNNNNGT");
    CHECK(reverse_complement_iupac("gaattc") == "GAATTC");  // either case in, upper case out
    CHECK(not reverse_complement_iupac("GAUC").has_value());
}
