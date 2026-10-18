#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_pair.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <charconv>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace {

using primerpicker::analysis::legacy_pair_values;
using primerpicker::analysis::PrimerPair;
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

// Reference values: PRIMER.exe, Analyze -> Primer Pair
// (tests/legacy/pairs/pair_screens.tsv, recorded 2026-10-18)

TEST_CASE("A:B dimers take the larger 3'-end count", "[analysis][legacy]") {
    auto const primer_a = nucleotides("GGAGGCGGAATGTAAACAGAAGG");
    auto const primer_b = nucleotides("ATGCTTAGGTGGATACCTTCT");
    auto const values = legacy_pair_values(PrimerPair{.first = primer_a, .second = primer_b});
    CHECK(values.dimer_3prime_matches == 6);  // PRIMER.exe
    CHECK(values.dimer_adjacent_bases == 6);  // PRIMER.exe
}

TEST_CASE("balance differences", "[analysis][legacy]") {
    auto const primer_a = nucleotides("CCCTGTCATACCAATCTACCCC");
    auto const primer_b = nucleotides("TGTTATGCGCGTGGGG");
    auto const values = legacy_pair_values(PrimerPair{.first = primer_a, .second = primer_b});
    CHECK(values.length_difference == 6);    // PRIMER.exe: "6 base difference"
    CHECK(values.gc_difference_percent == 8);  // PRIMER.exe: 54 and 62
    CHECK(values.tm_difference_celsius == 3);  // PRIMER.exe: 73 and 70
    CHECK(values.dimer_3prime_matches == 4);   // PRIMER.exe
    CHECK(values.dimer_adjacent_bases == 4);   // PRIMER.exe
}

TEST_CASE("pair values work at compile time", "[analysis][legacy][constexpr]") {
    // ATAATATTATATAGTGTTTT and TTTTATGGAAATAAAACAAA: A:B 3' end 6, any 6 (PRIMER.exe)
    constexpr auto values = [] {
        auto const primer_a = to_nucleotides("ATAATATTATATAGTGTTTT").value();
        auto const primer_b = to_nucleotides("TTTTATGGAAATAAAACAAA").value();
        auto const pair = legacy_pair_values(PrimerPair{.first = primer_a, .second = primer_b});
        return (pair.dimer_3prime_matches * 100) + pair.dimer_adjacent_bases;
    }();
    STATIC_CHECK(values == 606);
}

namespace {

struct Recording {
    std::string primer_a;
    std::string primer_b;
    std::vector<std::string> fields;
};

auto read_recordings(std::filesystem::path const& path) -> std::vector<Recording> {
    constexpr auto column_count = 26U;
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
        recordings.push_back({.primer_a = fields[0], .primer_b = fields[1], .fields = fields});
    }
    return recordings;
}

// a recorded value without its "!" flag
auto value(std::string text) -> int {
    if (text.ends_with('!')) {
        text.pop_back();
    }
    auto number = 0;
    auto const* const last = std::next(text.data(), std::ssize(text));
    auto const [end, error] = std::from_chars(text.data(), last, number);
    if (error != std::errc{} or end != last) {
        throw std::runtime_error{"unexpected number: " + text};
    }
    return number;
}

struct Columns {
    static constexpr auto length_a = 2U;
    static constexpr auto length_b = 3U;
    static constexpr auto gc_a = 5U;
    static constexpr auto gc_b = 6U;
    static constexpr auto tm_a = 8U;
    static constexpr auto tm_b = 9U;
    static constexpr auto dimer_3prime_ab = 22U;
    static constexpr auto dimer_adjacent_ab = 25U;
};

auto difference(Recording const& recording, std::size_t const column_a, std::size_t const column_b) -> int {
    auto const lhs = value(recording.fields.at(column_a));
    auto const rhs = value(recording.fields.at(column_b));
    return lhs > rhs ? lhs - rhs : rhs - lhs;
}

// the pair values, as integers that compare with ==
struct Values {
    int length_difference;
    int gc_difference_percent;
    int tm_difference_celsius;
    int dimer_3prime_matches;
    int dimer_adjacent_bases;

    friend auto operator==(Values const&, Values const&) -> bool = default;
};

auto recorded(Recording const& recording) -> Values {
    return {.length_difference = difference(recording, Columns::length_a, Columns::length_b),
            .gc_difference_percent = difference(recording, Columns::gc_a, Columns::gc_b),
            .tm_difference_celsius = difference(recording, Columns::tm_a, Columns::tm_b),
            .dimer_3prime_matches = value(recording.fields.at(Columns::dimer_3prime_ab)),
            .dimer_adjacent_bases = value(recording.fields.at(Columns::dimer_adjacent_ab))};
}

auto computed(Recording const& recording) -> Values {
    auto const primer_a = nucleotides(recording.primer_a);
    auto const primer_b = nucleotides(recording.primer_b);
    auto const values = legacy_pair_values(PrimerPair{.first = primer_a, .second = primer_b});
    return {.length_difference = values.length_difference,
            .gc_difference_percent = values.gc_difference_percent,
            .tm_difference_celsius = values.tm_difference_celsius,
            .dimer_3prime_matches = static_cast<int>(values.dimer_3prime_matches),
            .dimer_adjacent_bases = static_cast<int>(values.dimer_adjacent_bases)};
}

}  // namespace

TEST_CASE("legacy pair values match every recorded pair", "[analysis][legacy]") {
    auto const recordings = read_recordings(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "pairs" / "pair_screens.tsv");
    REQUIRE(recordings.size() >= 36);
    for (auto const& recording : recordings) {
        INFO("primers " << recording.primer_a << " " << recording.primer_b);
        CHECK(computed(recording) == recorded(recording));
    }
}
