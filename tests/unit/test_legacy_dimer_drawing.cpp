#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_dimer_drawing.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using primerpicker::analysis::legacy_dimer_drawing;
using primerpicker::analysis::PrimerPair;
using primerpicker::sequence::Nucleotide;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

struct RecordedDuplex {
    std::string first;
    std::string second;
    std::string top;
    std::string marks;
    std::string bottom;
};

// multi/multi_primer_duplexes.tsv: columns selection, list, row, column,
// duplex_5prime, duplex_pairs, duplex_3prime; the list holds the
// primers of fields A to D in order (selection 1234)
auto recorded_duplexes() -> std::vector<RecordedDuplex> {
    constexpr auto column_count = 7U;
    auto duplexes = std::vector<RecordedDuplex>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "multi" / "multi_primer_duplexes.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        if (fields.size() != column_count or fields[0] != "1234") {
            throw std::runtime_error{"unexpected line: " + line};
        }
        auto listed = std::vector<std::string>{};
        auto list_stream = std::istringstream{fields[1]};
        for (auto primer = std::string{}; std::getline(list_stream, primer, ',');) {
            listed.push_back(primer);
        }
        auto const field = [&listed](std::string const& name) -> std::string const& {
            return listed.at(static_cast<std::size_t>(name.at(0) - 'A'));
        };
        duplexes.push_back({.first = field(fields[2]),
                            .second = field(fields[3]),
                            .top = fields[4],
                            .marks = fields[5],
                            .bottom = fields[6]});
    }
    return duplexes;
}

void check_duplex(RecordedDuplex const& recorded) {
    CAPTURE(recorded.first, recorded.second);
    auto const first = nucleotides(recorded.first);
    auto const second = nucleotides(recorded.second);
    auto const drawing = legacy_dimer_drawing(PrimerPair{.first = first, .second = second});
    CHECK(drawing.top == recorded.top);
    CHECK(drawing.marks == recorded.marks);
    CHECK(drawing.bottom == recorded.bottom);
}

}  // namespace

// Reference values: PRIMER.exe, the dimer page of the Multi Primer
// Analysis (tests/legacy/multi/multi_primer_duplexes.tsv, recorded
// 2026-10-19)

TEST_CASE("dimer drawings match every recorded duplex", "[analysis][legacy]") {
    auto const duplexes = recorded_duplexes();
    REQUIRE(duplexes.size() == 32);
    for (auto const& recorded : duplexes) {
        check_duplex(recorded);
    }
}

TEST_CASE("the second primer's 3' end is drawn when its count is larger", "[analysis][legacy]") {
    // PRIMER.exe (multi_primer_duplexes.tsv): B GATCCATGCAATGCTAGCTA and
    // C ACGTTGCAGGTCAAGCTTGC, "3' end" 4: C's four 3'-most bases pair
    // with B
    check_duplex({.first = "GATCCATGCAATGCTAGCTA",
                  .second = "ACGTTGCAGGTCAAGCTTGC",
                  .top = "        5'-GATCCATGCAATGCTAGCTA-3'",
                  .marks = "                  ::::",
                  .bottom = "               3'-CGTTCGAACTGGACGTTGCA-5'"});
}

TEST_CASE("among equal counts, the largest partner position is drawn", "[analysis][legacy]") {
    // PRIMER.exe (multi_primer_duplexes.tsv): A CTCAAAGGCTCAGACATCAG on
    // B GATCCATGCAATGCTAGCTA: two pairs facing B's position 13 or 17; the
    // drawing shows 17
    check_duplex({.first = "CTCAAAGGCTCAGACATCAG",
                  .second = "GATCCATGCAATGCTAGCTA",
                  .top = "        5'-CTCAAAGGCTCAGACATCAG-3'",
                  .marks = "                             ::",
                  .bottom = "                         3'-ATCGATCGTAACGTACCTAG-5'"});
}
