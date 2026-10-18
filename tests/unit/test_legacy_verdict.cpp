#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using primerpicker::analysis::evaluate;
using primerpicker::analysis::legacy_gc_percent;
using primerpicker::analysis::LegacyActualValues;
using primerpicker::analysis::LegacyCriteria;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::to_nucleotides;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

}  // namespace

// Reference values: PRIMER.exe, built-in defaults
// (tests/legacy/README.md, "Built-in defaults" and "Primer criteria")

TEST_CASE("default criteria are those of PRIMER.exe", "[analysis][legacy]") {
    constexpr auto criteria = LegacyCriteria{};
    STATIC_CHECK(criteria.gc_min_percent == 50);                // PRIMER.exe
    STATIC_CHECK(criteria.gc_max_percent == 60);                // PRIMER.exe
    STATIC_CHECK(criteria.tm_min_celsius == 55);                // PRIMER.exe
    STATIC_CHECK(criteria.tm_max_celsius == 80);                // PRIMER.exe
    STATIC_CHECK(criteria.hairpin_cutoff.kcal_per_mol >= 0.0);  // PRIMER.exe: 0.0 kcal
    STATIC_CHECK(criteria.hairpin_cutoff.kcal_per_mol <= 0.0);
    STATIC_CHECK(criteria.reject_dimer_3prime_matches == 3);    // PRIMER.exe
    STATIC_CHECK(criteria.reject_dimer_adjacent_bases == 7);    // PRIMER.exe
    STATIC_CHECK(criteria.reject_run_length == 3);              // PRIMER.exe
    STATIC_CHECK(criteria.reject_gc_run_3prime == 3);           // PRIMER.exe
}

TEST_CASE("GC percent is truncated", "[analysis][legacy]") {
    CHECK(legacy_gc_percent(nucleotides("GTAAAACGACGGCCAGT")) == 52);  // PRIMER.exe: 9 of 17
    CHECK(legacy_gc_percent(nucleotides("AAGATTTA")) == 12);           // PRIMER.exe: 1 of 8
    CHECK(legacy_gc_percent(nucleotides("CCCGTCCC")) == 87);           // PRIMER.exe: 7 of 8
}

TEST_CASE("verdict of the M13 primer", "[analysis][legacy]") {
    // GTAAAACGACGGCCAGT, Inspect screen: all criteria met but runs
    constexpr auto actual = LegacyActualValues{.gc_percent = 52,
                                               .tm_celsius = 64,
                                               .hairpin = std::nullopt,
                                               .dimer_3prime_matches = 2,
                                               .dimer_adjacent_bases = 4,
                                               .longest_run = 4,
                                               .gc_run_3prime = 0};
    constexpr auto verdict = evaluate(actual, LegacyCriteria{});
    STATIC_CHECK(verdict.gc);                // PRIMER.exe: YES
    STATIC_CHECK(verdict.tm);                // PRIMER.exe: YES
    STATIC_CHECK(verdict.no_hairpin);        // PRIMER.exe: YES
    STATIC_CHECK(verdict.no_3prime_dimer);   // PRIMER.exe: YES
    STATIC_CHECK(verdict.no_dimer);          // PRIMER.exe: YES
    STATIC_CHECK(not verdict.no_runs);       // PRIMER.exe: NO
    STATIC_CHECK(verdict.no_3prime_gc_run);  // PRIMER.exe: YES
    STATIC_CHECK(not verdict.all());         // PRIMER.exe: MEETS CRITERIA SELECTED: NO
}

TEST_CASE("verdict of a primer with a hairpin and dimers", "[analysis][legacy]") {
    // AAAACCCCGGGGTTTTA, Inspect screen
    constexpr auto actual = LegacyActualValues{
        .gc_percent = 47,
        .tm_celsius = 68,
        .hairpin = primerpicker::thermo::KilocaloriesPerMole{-1.8},
        .dimer_3prime_matches = 2,
        .dimer_adjacent_bases = 16,
        .longest_run = 4,
        .gc_run_3prime = 0};
    constexpr auto verdict = evaluate(actual, LegacyCriteria{});
    STATIC_CHECK(not verdict.gc);          // PRIMER.exe: NO (47 < 50)
    STATIC_CHECK(verdict.tm);              // PRIMER.exe: YES
    STATIC_CHECK(not verdict.no_hairpin);  // PRIMER.exe: NO
    STATIC_CHECK(verdict.no_3prime_dimer); // PRIMER.exe: YES
    STATIC_CHECK(not verdict.no_dimer);    // PRIMER.exe: NO
    STATIC_CHECK(not verdict.no_runs);     // PRIMER.exe: NO
    STATIC_CHECK(verdict.no_3prime_gc_run);// PRIMER.exe: YES
}

TEST_CASE("range limits are inclusive", "[analysis][legacy]") {
    // PRIMER.exe shows YES for GC 50 and 60 %, and for Tm 55 and 80 °C
    // (recorded Inspect screens)
    auto actual = LegacyActualValues{.gc_percent = 50,
                                     .tm_celsius = 55,
                                     .hairpin = std::nullopt,
                                     .dimer_3prime_matches = 0,
                                     .dimer_adjacent_bases = 0,
                                     .longest_run = 0,
                                     .gc_run_3prime = 0};
    CHECK(evaluate(actual, LegacyCriteria{}).all());
    actual.gc_percent = 60;
    actual.tm_celsius = 80;
    CHECK(evaluate(actual, LegacyCriteria{}).all());
    actual.longest_run = 3;  // "reject runs of >= 3"
    CHECK_FALSE(evaluate(actual, LegacyCriteria{}).all());
}

namespace {

using primerpicker::analysis::legacy_actual_values;

}  // namespace

TEST_CASE("actual values of the M13 primer", "[analysis][legacy]") {
    // GTAAAACGACGGCCAGT, Inspect screen (PRIMER.exe)
    auto const actual = legacy_actual_values(nucleotides("GTAAAACGACGGCCAGT"));
    CHECK(actual.gc_percent == 52);
    CHECK(actual.tm_celsius == 64);
    CHECK_FALSE(actual.hairpin.has_value());  // "--"
    CHECK(actual.dimer_3prime_matches == 2);
    CHECK(actual.dimer_adjacent_bases == 4);
    CHECK(actual.longest_run == 4);
    CHECK(actual.gc_run_3prime == 0);
    CHECK_FALSE(evaluate(actual, LegacyCriteria{}).all());  // MEETS CRITERIA SELECTED: NO
}

TEST_CASE("the Tm is 0 when it is not positive", "[analysis][legacy]") {
    // AAGATTTA: "?" on the Direct Entry line, 0 on the Inspect screen (PRIMER.exe)
    CHECK(legacy_actual_values(nucleotides("AAGATTTA")).tm_celsius == 0);
}

TEST_CASE("actual values work at compile time", "[analysis][legacy][constexpr]") {
    // GGCCTCCGCC (PRIMER.exe, inspect.tsv): hairpin -0.4 kcal, 3' GC run 5
    constexpr auto gc_run = [] { return legacy_actual_values(to_nucleotides("GGCCTCCGCC").value()).gc_run_3prime; }();
    STATIC_CHECK(gc_run == 5);
}
