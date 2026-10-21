#include "core/analysis/legacy_hairpin.hpp"
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

using primerpicker::analysis::legacy_hairpin;
using primerpicker::analysis::legacy_hairpin_drawing;
using primerpicker::sequence::Nucleotide;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

struct RecordedDrawing {
    std::string primer;
    std::string top;     // hairpin_5prime
    std::string marks;   // hairpin_pairs
    std::string bottom;  // hairpin_3prime
};

// columns: primer, dimers_3prime, dimers_any, runs, gc_run_3prime,
// meets_no_hairpins, hairpin_kcal, hairpin_5prime, hairpin_pairs,
// hairpin_3prime (the drawing lines as on the screen, leading spaces
// included; "-" when there is no hairpin)
auto recorded_drawings() -> std::vector<RecordedDrawing> {
    constexpr auto column_count = 10U;
    auto drawings = std::vector<RecordedDrawing>{};
    auto const directory = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "inspect";
    // and the tandem G.T stems of 2026-10-21 (approved the same day)
    for (auto const* const name : {"inspect.tsv", "hairpin_design.tsv", "inspect_2.tsv",
                                   "../hairpin_tandem_gt_2/hairpin_tandem_gt_2.tsv"}) {
        auto file = std::ifstream{directory / name};
        for (auto line = std::string{}; std::getline(file, line);) {
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
            if (fields[6] != "-") {
                drawings.push_back({.primer = fields[0], .top = fields[7], .marks = fields[8], .bottom = fields[9]});
            }
        }
    }
    return drawings;
}

// the drawing of TEXT against the three recorded lines
void check_drawing(std::string const& text, RecordedDrawing const& expected) {
    CAPTURE(text);
    auto const primer = nucleotides(text);
    auto const hairpin = legacy_hairpin(primer);
    REQUIRE(hairpin.has_value());
    auto const drawing = legacy_hairpin_drawing(primer, hairpin.value_or(primerpicker::analysis::Hairpin{}));
    CHECK(drawing.top == expected.top);
    CHECK(drawing.marks == expected.marks);
    CHECK(drawing.bottom == expected.bottom);
}

auto shown_energy(std::string const& text) -> double {
    auto stream = std::istringstream{text};
    auto value = 0.0;
    stream >> value;
    if (stream.fail() or not stream.eof()) {
        throw std::runtime_error{"unexpected energy: " + text};
    }
    return value;
}

auto computed_energy(std::string const& text) -> double {
    auto const hairpin = legacy_hairpin(nucleotides(text));
    if (not hairpin.has_value()) {
        throw std::runtime_error{"no hairpin: " + text};
    }
    return primerpicker::analysis::legacy_displayed_energy(*hairpin).kcal_per_mol;
}

// multi/multi_primer_hairpins.tsv: columns selection, list, field,
// hairpin_kcal, hairpin_5prime, hairpin_pairs, hairpin_3prime; the list
// holds the primers of fields A to D in order (selection 1234)
struct PageRow {
    std::string primer;
    std::string energy;
    RecordedDrawing drawing;
};

auto hairpin_page_rows() -> std::vector<PageRow> {
    constexpr auto column_count = 7U;
    auto rows = std::vector<PageRow>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "multi" / "multi_primer_hairpins.tsv"};
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
        auto const& primer = listed.at(static_cast<std::size_t>(fields[2].at(0) - 'A'));
        rows.push_back({.primer = primer,
                        .energy = fields[3],
                        .drawing = {.primer = primer, .top = fields[4], .marks = fields[5], .bottom = fields[6]}});
    }
    return rows;
}

}  // namespace

// Reference values: PRIMER.exe, the hairpin drawing of the Primer
// Analysis (Inspect) screen (tests/legacy/inspect/*.tsv, recorded
// 2026-10-17 and 2026-10-18; tests/legacy/hairpin_tandem_gt_2/, recorded
// 2026-10-21)

TEST_CASE("hairpin drawings match every recorded drawing", "[analysis][legacy]") {
    auto const drawings = recorded_drawings();
    REQUIRE(drawings.size() == 285);
    for (auto const& recorded : drawings) {
        check_drawing(recorded.primer, recorded);
    }
}

TEST_CASE("a hairpin drawing: stem, then half the loop on each line", "[analysis][legacy]") {
    // PRIMER.exe (inspect/hairpin_design.tsv): GACTTTTTAGTC, stem of four
    // pairs, loop of four bases
    auto const primer = nucleotides("GACTTTTTAGTC");
    auto const hairpin = legacy_hairpin(primer);
    REQUIRE(hairpin.has_value());
    auto const drawing = legacy_hairpin_drawing(primer, hairpin.value_or(primerpicker::analysis::Hairpin{}));
    CHECK(drawing.top == "             5'-GACTTT");
    CHECK(drawing.marks == "                ::::  )");
    CHECK(drawing.bottom == "             3'-CTGATT");
}

TEST_CASE("nucleotides convert back to letters", "[sequence]") {
    CHECK(primerpicker::sequence::to_letter(Nucleotide::adenine) == 'A');
    CHECK(primerpicker::sequence::to_letter(Nucleotide::cytosine) == 'C');
    CHECK(primerpicker::sequence::to_letter(Nucleotide::guanine) == 'G');
    CHECK(primerpicker::sequence::to_letter(Nucleotide::thymine) == 'T');
}

TEST_CASE("the middle base of an odd loop is on both lines", "[analysis][legacy]") {
    // PRIMER.exe (inspect/*.tsv): ACTAGTAGCCCACCTACTGTAT, loop of five
    // bases: C (position 11) ends the 5' line and starts the 3' line
    auto const primer = nucleotides("ACTAGTAGCCCACCTACTGTAT");
    auto const hairpin = legacy_hairpin(primer);
    REQUIRE(hairpin.has_value());
    auto const drawing = legacy_hairpin_drawing(primer, hairpin.value_or(primerpicker::analysis::Hairpin{}));
    CHECK(drawing.top == "             5'-ACTAGTAGCCC");
    CHECK(drawing.bottom == "            3'-TATGTCATCCAC");
}

TEST_CASE("a long 3' part is cut at the left edge", "[analysis][legacy]") {
    // PRIMER.exe (inspect/*.tsv): the 3' line loses its label and its
    // 3'-most base, which would be left of column 1
    auto const primer = nucleotides("GCACCTGCCAAAGGCGGCGGGGTAAATACCGCCCCGGCGC");
    auto const hairpin = legacy_hairpin(primer);
    REQUIRE(hairpin.has_value());
    auto const drawing = legacy_hairpin_drawing(primer, hairpin.value_or(primerpicker::analysis::Hairpin{}));
    CHECK(drawing.bottom == " GCGGCCCCGCCATAAATGGGGCGGCGG");
}

TEST_CASE("the Multi Primer hairpin page draws the same hairpins", "[analysis][legacy]") {
    // PRIMER.exe (multi/multi_primer_hairpins.tsv, recorded 2026-10-19)
    auto const rows = hairpin_page_rows();
    REQUIRE(rows.size() == 4);
    constexpr auto tenth = 0.05;
    for (auto const& row : rows) {
        check_drawing(row.primer, row.drawing);
        auto const difference = computed_energy(row.primer) - shown_energy(row.energy);
        CHECK((difference < tenth and difference > -tenth));
    }
}
