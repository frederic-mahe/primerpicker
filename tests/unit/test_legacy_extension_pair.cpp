#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_pair.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <optional>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

using primerpicker::analysis::LegacyCriteria;
using primerpicker::analysis::LegacyMultiPrimerFields;
using primerpicker::analysis::LegacyPrimer;
using primerpicker::analysis::LegacyPrimerPair;
using primerpicker::sequence::Nucleotide;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence"};
    }
    return *converted;
}

// the value; throws (and fails the test) when there is none
template <typename Value>
auto checked(std::optional<Value> const& optional) -> Value {
    if (not optional.has_value()) {
        throw std::runtime_error{"no value"};
    }
    return *optional;
}

}  // namespace

// Reference: PRIMER.exe, Analyze -> Primer Pair for the primer made at 51
// of TEST600 with the 5' extension GAATTC (20 bases in all) and the
// primer at 400(C) (tests/legacy/extension/ext_pair.screens.txt,
// recorded 2026-10-22): the values of the whole primer, the Tm without
// the extension
TEST_CASE("Primer Pair with a 5' extension", "[analysis][legacy]") {
    auto const primer_a = nucleotides("GAATTCCTCAAAGGCTCAGA");
    auto const primer_b = nucleotides("TGTATCCAGGTAATGCAAGG");
    auto const extended_a = LegacyPrimer{.sequence = primer_a, .extension_length = 6};
    auto const plain_b = LegacyPrimer{.sequence = primer_b, .extension_length = 0};
    auto const actual_a = primerpicker::analysis::legacy_actual_values(extended_a);
    CHECK(actual_a.gc_percent == 45);
    CHECK(actual_a.tm_celsius == 47);
    auto const pair = primerpicker::analysis::legacy_pair_values(LegacyPrimerPair{.first = extended_a, .second = plain_b});
    CHECK(pair.length_difference == 0);      // 20 and 20
    CHECK(pair.tm_difference_celsius == 19);  // "19°C difference"
    CHECK(pair.gc_difference_percent == 0);
    CHECK(pair.dimer_3prime_matches == 3);  // A:B "3!"
    CHECK(pair.dimer_adjacent_bases == 3);
}

// Reference: PRIMER.exe, Analyze -> Multi Primer for M13EXT.PD2 (GAATTC +
// GTAAAACGACGGCCAGT) and AACAGCTATGACCATG
// (tests/legacy/extension/ext_multi.screens.txt, recorded 2026-10-22)
TEST_CASE("Multi Primer with a 5' extension", "[analysis][legacy]") {
    auto const primer_a = nucleotides("GAATTCGTAAAACGACGGCCAGT");
    auto const primer_b = nucleotides("AACAGCTATGACCATG");
    auto const fields =
        LegacyMultiPrimerFields{LegacyPrimer{.sequence = primer_a, .extension_length = 6},
                                LegacyPrimer{.sequence = primer_b, .extension_length = 0}, std::nullopt, std::nullopt};
    auto const analysis = primerpicker::analysis::legacy_multi_primer(fields, LegacyCriteria{});
    auto const column_a = checked(analysis.columns.at(0));
    CHECK(column_a.length == 23);
    CHECK(column_a.actual.gc_percent == 47);
    CHECK(column_a.cautions.gc);
    CHECK(column_a.actual.tm_celsius == 64);
    CHECK_FALSE(column_a.cautions.tm);
    CHECK(column_a.actual.gc_run_3prime == 0);
    CHECK(column_a.actual.longest_run == 4);
    CHECK(column_a.cautions.runs);
    CHECK_FALSE(column_a.cautions.hairpin);
    auto const self_a = checked(analysis.dimers.at(0).at(0));
    CHECK(self_a.dimer_3prime_matches == 2);
    CHECK(self_a.dimer_adjacent_bases == 6);
    auto const a_b = checked(analysis.dimers.at(0).at(1));
    CHECK(a_b.dimer_3prime_matches == 2);
    CHECK(a_b.dimer_adjacent_bases == 2);
}
