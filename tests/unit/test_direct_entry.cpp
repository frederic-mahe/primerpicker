#include "core/sequence/legacy_direct_entry.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string_view>

using primerpicker::sequence::direct_entry_inserts;
using primerpicker::sequence::direct_entry_max_length;

// Reference: PRIMER.exe, Primer/Create -> Direct Entry
// (tests/legacy/README.md, "Primer Tm", recorded 2026-10-17): "The
// primer field holds 50 bases; further keys are ignored."
TEST_CASE("the Direct Entry field holds 50 bases", "[direct entry]") {
    STATIC_REQUIRE(direct_entry_max_length == 50);
}

// Same reference: "Letters other than A, C, G and T (N, R, U) are not
// inserted"
TEST_CASE("Direct Entry inserts A, C, G and T only", "[direct entry]") {
    for (auto const letter : std::string_view{"ACGT"}) {
        CHECK(direct_entry_inserts(letter));
    }
    for (auto const letter : std::string_view{"NRU"}) {
        CHECK_FALSE(direct_entry_inserts(letter));
    }
}
