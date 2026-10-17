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

auto to_int(std::string const& text) -> std::optional<int> {
    auto value = 0;
    auto const* const last = std::next(text.data(), std::ssize(text));
    auto const [end, error] = std::from_chars(text.data(), last, value);
    if (error != std::errc{} or end != last) {
        return std::nullopt;
    }
    return value;
}

struct Recording {
    std::string primer;
    std::string shown_tm;
};

// columns: typed_extension, typed_primer, shown_length,
// shown_gc_percent, gc_flag, shown_tm_celsius, tm_flag
auto read_recordings(std::filesystem::path const& path) -> std::vector<Recording> {
    constexpr auto column_count = 7U;
    constexpr auto primer_column = 1U;
    constexpr auto tm_column = 5U;
    auto recordings = std::vector<Recording>{};
    auto file = std::ifstream{path};
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
        assert(fields.size() == column_count);
        recordings.push_back({fields[primer_column], fields[tm_column]});
    }
    return recordings;
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

namespace {

// Entries that the program edited are skipped: longer than its 50-base
// field, or with letters it rejects
auto is_unedited(std::string const& primer) -> bool {
    constexpr auto field_length = 50U;
    return to_nucleotides(primer).has_value() and primer.size() <= field_length;
}

auto computed_tm(std::string const& primer) -> std::optional<int> {
    return legacy_displayed_tm(legacy_primer_tm(nucleotides(primer)));
}

auto shown_tm(std::string const& shown) -> std::optional<int> {
    auto const value = shown == "?" ? undefined : to_int(shown);
    assert(shown == "?" or value.has_value());
    return value;
}

}  // namespace

TEST_CASE("legacy Tm matches every recorded primer", "[thermo][legacy]") {
    auto const recordings =
        read_recordings(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "tm" / "direct_entry_tm.tsv");
    REQUIRE(recordings.size() == 186);
    auto checked = 0;
    for (auto const& [primer, shown] : recordings) {
        if (not is_unedited(primer)) {
            continue;
        }
        INFO("primer " << primer << ", shown " << shown);
        CHECK(computed_tm(primer) == shown_tm(shown));
        ++checked;
    }
    CHECK(checked == 180);
}
