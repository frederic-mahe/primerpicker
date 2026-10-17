#include "core/sequence/nucleotide.hpp"
#include "core/thermo/legacy_tm.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cassert>
#include <charconv>
#include <iterator>
#include <system_error>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace {

using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::to_nucleotides;
using primerpicker::thermo::legacy_displayed_tm;
using primerpicker::thermo::legacy_primer_tm;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = to_nucleotides(letters);
    assert(converted.has_value());
    return *converted;
}

auto displayed_tm(std::string_view const primer) -> std::optional<int> {
    return legacy_displayed_tm(legacy_primer_tm(nucleotides(primer)));
}

constexpr auto undefined = std::optional<int>{};

auto to_int(std::string_view const text) -> std::optional<int> {
    auto value = 0;
    auto const [end, error] = std::from_chars(text.data(), std::next(text.data(), std::ssize(text)), value);
    if (error != std::errc{} or end != std::next(text.data(), std::ssize(text))) {
        return std::nullopt;
    }
    return value;
}

}  // namespace

// Reference values: PRIMER.exe, Primer/Create -> Direct Entry, with
// built-in defaults (tests/legacy/tm/direct_entry_tm.tsv); "?" is
// shown when the Tm is not positive.

TEST_CASE("legacy Tm of homopolymers and repeats", "[thermo][legacy]") {
    CHECK(displayed_tm("AAAAAAAAAAAA") == 35);                    // PRIMER.exe
    CHECK(displayed_tm("ATATATATATATATATATATATATATATAT") == 51);  // PRIMER.exe
    CHECK(displayed_tm("CGCGCGCGCGCGCGCGCGCGCGCGCGCGCG") == 119); // PRIMER.exe
}

TEST_CASE("legacy Tm of mixed primers", "[thermo][legacy]") {
    CHECK(displayed_tm("CCCGTCCC") == 28);                   // PRIMER.exe
    CHECK(displayed_tm("AGCCCGTAACGATACAT") == 60);          // PRIMER.exe
    CHECK(displayed_tm("TCGGGGGGTTATTATAGGTGAAGTT") == 77);  // PRIMER.exe
    CHECK(displayed_tm("CCAGTTAACATACAGGATTACTAACCGTTGCGGCATAAGTTGGGTTGAAT") == 96);  // PRIMER.exe
}

TEST_CASE("legacy Tm is undefined when not positive", "[thermo][legacy]") {
    CHECK(displayed_tm("AAGATTTA") == undefined);   // PRIMER.exe: Tm about -6.9
    CHECK(displayed_tm("ACTATTGAA") == undefined);  // PRIMER.exe: Tm between 0 and 1
}

TEST_CASE("legacy Tm matches every recorded primer", "[thermo][legacy]") {
    // columns: typed_extension, typed_primer, shown_length,
    // shown_gc_percent, gc_flag, shown_tm_celsius, tm_flag
    auto const path = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "tm" / "direct_entry_tm.tsv";
    auto file = std::ifstream{path};
    REQUIRE(file.is_open());
    auto checked = 0;
    auto line = std::string{};
    while (std::getline(file, line)) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        REQUIRE(fields.size() == 7);
        auto const& primer = fields[1];
        auto const& shown = fields[5];
        auto const converted = to_nucleotides(primer);
        // skip entries that the program edited: longer than its
        // 50-base field, or with letters it rejects
        constexpr auto field_length = 50U;
        if (not converted.has_value() or primer.size() > field_length) {
            continue;
        }
        auto const expected = shown == "?" ? undefined : to_int(shown);
        REQUIRE((shown == "?" or expected.has_value()));
        INFO("primer " << primer << ", shown " << shown);
        CHECK(legacy_displayed_tm(legacy_primer_tm(*converted)) == expected);
        ++checked;
    }
    CHECK(checked == 180);
}
