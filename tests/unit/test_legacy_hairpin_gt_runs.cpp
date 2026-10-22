#include "core/analysis/legacy_hairpin.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
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

// hairpin_gt_runs/hairpin_gt_runs.tsv: columns primer, dimers_3prime,
// dimers_any, runs, gc_run_3prime, meets_no_hairpins, hairpin_kcal,
// hairpin_5prime, hairpin_pairs, hairpin_3prime
auto recorded_rows() -> std::vector<std::vector<std::string>> {
    constexpr auto column_count = 10U;
    auto rows = std::vector<std::vector<std::string>>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "hairpin_gt_runs" / "hairpin_gt_runs.tsv"};
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
        rows.push_back(std::move(fields));
    }
    return rows;
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

// the energy and drawing of a recorded row against primerpicker's
void check_row(std::vector<std::string> const& row) {
    CAPTURE(row[0]);
    auto const primer = nucleotides(row[0]);
    auto const hairpin = legacy_hairpin(primer);
    REQUIRE(hairpin.has_value());
    auto const found = hairpin.value_or(primerpicker::analysis::Hairpin{});
    // tenths of kcal: an explicit tolerance below the shown precision
    constexpr auto tolerance = 1e-9;
    auto const difference = primerpicker::analysis::legacy_displayed_energy(found).kcal_per_mol - shown_energy(row[6]);
    CHECK((difference < tolerance and difference > -tolerance));
    auto const drawing = legacy_hairpin_drawing(primer, found);
    CHECK(drawing.top == row[7]);
    CHECK(drawing.marks == row[8]);
    CHECK(drawing.bottom == row[9]);
}

}  // namespace

// Reference values: PRIMER.exe, the Hairpins screen of 32 primers whose
// stem holds a run of 3 or 4 G.T pairs (tests/legacy/hairpin_gt_runs/,
// recorded 2026-10-22)
TEST_CASE("hairpins with runs of G.T pairs: every recorded energy and drawing", "[analysis][legacy]") {
    auto const rows = recorded_rows();
    REQUIRE(rows.size() == 32);
    for (auto const& row : rows) {
        check_row(row);
    }
}
