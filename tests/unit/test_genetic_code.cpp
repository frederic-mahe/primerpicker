#include "core/sequence/genetic_code.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstddef>
#include <string_view>

namespace {

using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::translate;

// NCBI, "The Genetic Codes", translation table 1 (the standard code):
// amino acids of the 64 codons, the three bases running through T, C,
// A, G in turn (first base slowest); '*' is a stop codon
constexpr auto ncbi_table_1 = std::string_view{"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"};
constexpr auto ncbi_order = std::array{Nucleotide::thymine, Nucleotide::cytosine, Nucleotide::adenine,
                                       Nucleotide::guanine};

}  // namespace

TEST_CASE("the standard genetic code (NCBI table 1), stops shown as '-'", "[sequence]") {
    auto index = std::size_t{0};
    for (auto const first : ncbi_order) {
        for (auto const second : ncbi_order) {
            for (auto const third : ncbi_order) {
                auto const expected = ncbi_table_1[index] == '*' ? '-' : ncbi_table_1[index];
                CAPTURE(index);
                CHECK(translate(first, second, third) == expected);
                ++index;
            }
        }
    }
}

TEST_CASE("codons as PRIMER.exe translates them", "[sequence][legacy]") {
    // PRIMER.exe, Mutagenize Primer (tests/legacy/mutagenize/, recorded
    // 2026-10-19): CTC L, ATC I, TAA -, GAG E
    CHECK(translate(Nucleotide::cytosine, Nucleotide::thymine, Nucleotide::cytosine) == 'L');
    CHECK(translate(Nucleotide::adenine, Nucleotide::thymine, Nucleotide::cytosine) == 'I');
    CHECK(translate(Nucleotide::thymine, Nucleotide::adenine, Nucleotide::adenine) == '-');
    CHECK(translate(Nucleotide::guanine, Nucleotide::adenine, Nucleotide::guanine) == 'E');
}

TEST_CASE("translation works at compile time", "[sequence][constexpr]") {
    STATIC_CHECK(translate(Nucleotide::adenine, Nucleotide::thymine, Nucleotide::guanine) == 'M');
}
