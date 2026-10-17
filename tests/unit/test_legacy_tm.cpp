#include "core/sequence/nucleotide.hpp"
#include "core/thermo/legacy_tm.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <charconv>
#include <iterator>
#include <system_error>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
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
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
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
        if (fields.size() != column_count) {
            throw std::runtime_error{"unexpected line: " + line};
        }
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
    if (shown != "?" and not value.has_value()) {
        throw std::runtime_error{"unexpected Tm: " + shown};
    }
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

// Reference values: PRIMER.exe, Analyze -> Primer Pair -> AnalProd,
// on generated molecules (tests/legacy/pairs/, recorded 2026-10-17)

namespace {

using primerpicker::thermo::legacy_annealing_tm;
using primerpicker::thermo::legacy_product_tm;
using primerpicker::thermo::LegacyPairTms;

auto displayed_product_tm(std::string_view const product) -> std::optional<int> {
    return legacy_displayed_tm(legacy_product_tm(nucleotides(product)));
}

}  // namespace

TEST_CASE("legacy product Tm of recorded products", "[thermo][legacy]") {
    // TGRAD, primers at 961 and 1031C: 71 bp
    CHECK(displayed_product_tm("ATAATATTATATAGTGTTTTTTTGTCGTCTATGTATGATATTTTTATTATATTTGTTTTATTTCC"
                               "ATAAAA") == 56);  // PRIMER.exe
    // TGRAD, primers at 322 and 380C: 59 bp
    CHECK(displayed_product_tm("CACCAGCTAGTTATGTTCGTTGTGTGGACAGTCCCAACATCTGCGACTCCACTTCCCTA") ==
          69);  // PRIMER.exe
    // TGRAD, primers at 148 and 229C: 82 bp
    CHECK(displayed_product_tm("CGTAGTTTACTTACAGAACATAGGCCTGTAGACGTTACAAATATCGATCTTTATTACTTGACTG"
                               "GGGCCAATGTTATCTCTT") == 67);  // PRIMER.exe
}

TEST_CASE("legacy annealing temperature uses the lower primer Tm", "[thermo][legacy]") {
    // TEST600, primers at 51 and 400C: Tm 66 and 66, product Tm 77
    CHECK(legacy_displayed_tm(legacy_annealing_tm(LegacyPairTms{
              .primer_a_celsius = 66, .primer_b_celsius = 66, .product_celsius = 77})) == 52);  // PRIMER.exe
    // TGRAD, primers at 43 and 1078C: Tm 44 and 92, product Tm 78
    CHECK(legacy_displayed_tm(legacy_annealing_tm(LegacyPairTms{
              .primer_a_celsius = 44, .primer_b_celsius = 92, .product_celsius = 78})) == 46);  // PRIMER.exe
    // TGRAD, primers at 195 and 969C: Tm 63 and 49, product Tm 79
    CHECK(legacy_displayed_tm(legacy_annealing_tm(LegacyPairTms{
              .primer_a_celsius = 63, .primer_b_celsius = 49, .product_celsius = 79})) == 48);  // PRIMER.exe
    // TGRAD, primers at 961 and 1031C: Tm 47 and 59, product Tm 56
    CHECK(legacy_displayed_tm(legacy_annealing_tm(LegacyPairTms{
              .primer_a_celsius = 47, .primer_b_celsius = 59, .product_celsius = 56})) == 31);  // PRIMER.exe
}

TEST_CASE("legacy Tm functions work at compile time", "[thermo][legacy][constexpr]") {
    constexpr auto primer = std::array{Nucleotide::cytosine, Nucleotide::cytosine, Nucleotide::cytosine,
                                       Nucleotide::guanine, Nucleotide::thymine, Nucleotide::cytosine,
                                       Nucleotide::cytosine, Nucleotide::cytosine};
    // CCCGTCCC: 28 °C shown by PRIMER.exe (direct_entry_tm.tsv)
    STATIC_CHECK(legacy_primer_tm(primer).degrees >= 28.0);
    STATIC_CHECK(legacy_primer_tm(primer).degrees < 29.0);
    // TEST600, primers at 51 and 400C: Ta 52 °C shown by PRIMER.exe (pairs.tsv)
    constexpr auto pair_tms = LegacyPairTms{.primer_a_celsius = 66, .primer_b_celsius = 66, .product_celsius = 77};
    STATIC_CHECK(legacy_annealing_tm(pair_tms).degrees >= 52.0);
    STATIC_CHECK(legacy_annealing_tm(pair_tms).degrees < 53.0);
    // TGRAD, primers at 322 and 380C: product Tm 69 °C shown by PRIMER.exe (pairs.tsv)
    constexpr auto product_tm = [] {
        auto const product =
            to_nucleotides("CACCAGCTAGTTATGTTCGTTGTGTGGACAGTCCCAACATCTGCGACTCCACTTCCCTA").value();
        return legacy_product_tm(product).degrees;
    }();
    STATIC_CHECK(product_tm >= 69.0);
    STATIC_CHECK(product_tm < 70.0);
}
