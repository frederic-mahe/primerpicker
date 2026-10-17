#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <vector>

using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::to_nucleotide;
using primerpicker::sequence::to_nucleotides;

TEST_CASE("each of A, C, G and T is a nucleotide", "[sequence]") {
    CHECK(to_nucleotide('A') == Nucleotide::adenine);
    CHECK(to_nucleotide('C') == Nucleotide::cytosine);
    CHECK(to_nucleotide('G') == Nucleotide::guanine);
    CHECK(to_nucleotide('T') == Nucleotide::thymine);
}

TEST_CASE("lowercase letters are nucleotides too", "[sequence]") {
    CHECK(to_nucleotide('a') == Nucleotide::adenine);
    CHECK(to_nucleotide('c') == Nucleotide::cytosine);
    CHECK(to_nucleotide('g') == Nucleotide::guanine);
    CHECK(to_nucleotide('t') == Nucleotide::thymine);
}

TEST_CASE("other characters are not nucleotides", "[sequence]") {
    CHECK_FALSE(to_nucleotide('N').has_value());
    CHECK_FALSE(to_nucleotide('U').has_value());
    CHECK_FALSE(to_nucleotide(' ').has_value());
    CHECK_FALSE(to_nucleotide('\0').has_value());
}

TEST_CASE("a string of nucleotides converts in order", "[sequence]") {
    auto const expected = std::vector{Nucleotide::guanine, Nucleotide::adenine,
                                      Nucleotide::thymine, Nucleotide::cytosine};
    CHECK(to_nucleotides("GaTc") == expected);
}

TEST_CASE("an empty string converts to an empty sequence", "[sequence]") {
    CHECK(to_nucleotides("") == std::vector<Nucleotide>{});
}

TEST_CASE("a string with a non-nucleotide does not convert", "[sequence]") {
    CHECK_FALSE(to_nucleotides("GTAAAACGNCGG").has_value());
}

TEST_CASE("Watson-Crick complements", "[sequence]") {
    using primerpicker::sequence::complement;
    STATIC_CHECK(complement(Nucleotide::adenine) == Nucleotide::thymine);
    STATIC_CHECK(complement(Nucleotide::cytosine) == Nucleotide::guanine);
    STATIC_CHECK(complement(Nucleotide::guanine) == Nucleotide::cytosine);
    STATIC_CHECK(complement(Nucleotide::thymine) == Nucleotide::adenine);
}

TEST_CASE("nucleotide conversions work at compile time", "[sequence][constexpr]") {
    STATIC_CHECK(to_nucleotide('g') == Nucleotide::guanine);
    STATIC_CHECK(not to_nucleotide('N').has_value());
    STATIC_CHECK(to_nucleotides("GAT").value().size() == 3);
}
