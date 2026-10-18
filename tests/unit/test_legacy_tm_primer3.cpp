#include "core/sequence/nucleotide.hpp"
#include "core/thermo/legacy_tm.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Catch::Matchers::WithinAbs;
using primerpicker::sequence::to_nucleotides;
using primerpicker::thermo::legacy_primer_tm;

struct ReferenceTm {
    std::string primer;
    double tm_celsius;
};

// tests/data/primer3_legacy_tm.tsv: primer, then the Tm printed by
// Primer3's oligotm
auto reference_values() -> std::vector<ReferenceTm> {
    auto values = std::vector<ReferenceTm>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_TEST_DATA_DIR} / "primer3_legacy_tm.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto stream = std::istringstream{line};
        auto primer = std::string{};
        auto tm_celsius = 0.0;
        if (not(stream >> primer >> tm_celsius)) {
            throw std::runtime_error{"cannot read line: " + line};
        }
        values.push_back({.primer = std::move(primer), .tm_celsius = tm_celsius});
    }
    return values;
}

}  // namespace

// Reference values: Primer3 2.6.1, oligotm -tp 0 -sc 0 -dv 0 -n 0
// -d 0.0496293659554817 -mv 1021.0244559195794 (Breslauer 1986
// parameters; DNA and salt concentrations chosen to match Primer
// Designer's constants, see tests/data/primer3_legacy_tm.sh), for the
// primers recorded from PRIMER.exe (tests/legacy/tm/direct_entry_tm.tsv)

TEST_CASE("Legacy primer Tm equals Primer3's Breslauer Tm with matching settings", "[thermo][legacy][primer3]") {
    auto const references = reference_values();
    REQUIRE(references.size() == 159);
    // oligotm prints 6 decimals
    constexpr auto tolerance_celsius = 1e-5;
    for (auto const& [primer, tm_celsius] : references) {
        INFO(primer);
        auto const sequence = to_nucleotides(primer);
        if (not sequence.has_value()) {
            FAIL("not a nucleotide sequence");
            continue;
        }
        CHECK_THAT(legacy_primer_tm(*sequence).degrees, WithinAbs(tm_celsius, tolerance_celsius));
    }
}
