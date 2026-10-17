#include "core/sequence/nucleotide.hpp"
#include "core/thermo/breslauer_1986.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using Catch::Matchers::WithinAbs;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::to_nucleotides;
using primerpicker::thermo::breslauer_1986_sums;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

// sums are given with one decimal in the source
constexpr auto tolerance = 0.005;

}  // namespace

// Reference values: Primer3 2.6.1, oligotm -tp 0 -sc 0 -mv 1000 -dv 0
// -n 0, run at two DNA concentrations (1e9 and 1e-3 nM) to solve for
// dH and dS; Primer3's constant initiation entropy (-10.8 cal/K/mol)
// is removed. Signs are those of dissociation (positive).

TEST_CASE("Breslauer 1986 sums of a homopolymer", "[thermo]") {
    auto const sums = breslauer_1986_sums(nucleotides("AAAAAAAAAAAAAAAA"));
    CHECK_THAT(sums.enthalpy.kcal_per_mol, WithinAbs(136.5, tolerance));  // Primer3 2.6.1
    CHECK_THAT(sums.entropy.cal_per_kelvin_per_mol, WithinAbs(360.0, tolerance));  // Primer3 2.6.1
}

TEST_CASE("Breslauer 1986 sums of a mixed primer", "[thermo]") {
    auto const sums = breslauer_1986_sums(nucleotides("GTAAAACGACGGCCAGT"));
    CHECK_THAT(sums.enthalpy.kcal_per_mol, WithinAbs(135.4, tolerance));  // Primer3 2.6.1
    CHECK_THAT(sums.entropy.cal_per_kelvin_per_mol, WithinAbs(340.8, tolerance));  // Primer3 2.6.1
}

TEST_CASE("Breslauer 1986 sums cover all 16 pairs", "[thermo]") {
    // the 16 dinucleotide repeats of length 16 together use every pair
    struct Expected {
        std::string_view primer;
        double enthalpy_kcal_per_mol;
        double entropy_cal_per_kelvin_per_mol;
    };
    // Primer3 2.6.1 (see above)
    auto const cases = std::vector<Expected>{
        {"ACACACACACACACAC", 92.6, 228.7},  {"AGAGAGAGAGAGAGAG", 101.6, 260.9},
        {"ATATATATATATATAT", 110.8, 309.5}, {"CACACACACACACACA", 91.9, 224.3},
        {"CCCCCCCCCCCCCCCC", 165.0, 399.0}, {"CGCGCGCGCGCGCGCG", 172.9, 409.3},
        {"CTCTCTCTCTCTCTCT", 101.6, 260.9}, {"GAGAGAGAGAGAGAGA", 99.4, 253.6},
        {"GCGCGCGCGCGCGCGC", 172.1, 408.2}, {"GGGGGGGGGGGGGGGG", 165.0, 399.0},
        {"GTGTGTGTGTGTGTGT", 92.6, 228.7},  {"TATATATATATATATA", 108.2, 302.5},
        {"TCTCTCTCTCTCTCTC", 99.4, 253.6},  {"TGTGTGTGTGTGTGTG", 91.9, 224.3},
        {"TTTTTTTTTTTTTTTT", 136.5, 360.0},
    };
    for (auto const& expected : cases) {
        INFO(expected.primer);
        auto const sums = breslauer_1986_sums(nucleotides(expected.primer));
        CHECK_THAT(sums.enthalpy.kcal_per_mol, WithinAbs(expected.enthalpy_kcal_per_mol, tolerance));
        CHECK_THAT(sums.entropy.cal_per_kelvin_per_mol,
                   WithinAbs(expected.entropy_cal_per_kelvin_per_mol, tolerance));
    }
}

TEST_CASE("Breslauer 1986 sums work at compile time", "[thermo][constexpr]") {
    using primerpicker::sequence::Nucleotide;
    constexpr auto poly_a = std::array{Nucleotide::adenine, Nucleotide::adenine, Nucleotide::adenine};
    // two AA pairs: 2 x 9.1 kcal/mol (Breslauer 1986, Table 2)
    STATIC_CHECK(breslauer_1986_sums(poly_a).enthalpy.kcal_per_mol > 18.19);
    STATIC_CHECK(breslauer_1986_sums(poly_a).enthalpy.kcal_per_mol < 18.21);
}
