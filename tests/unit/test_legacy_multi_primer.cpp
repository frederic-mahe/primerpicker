#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using primerpicker::analysis::LegacyCriteria;
using primerpicker::analysis::legacy_multi_primer;
using primerpicker::analysis::LegacyMultiPrimer;
using primerpicker::analysis::LegacyMultiPrimerColumn;
using primerpicker::analysis::MultiPrimerFields;
using primerpicker::sequence::Nucleotide;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

auto read_table(std::string_view const name) -> std::vector<std::vector<std::string>> {
    auto rows = std::vector<std::vector<std::string>>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "multi" / name};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        rows.push_back(std::move(fields));
    }
    if (rows.empty()) {
        throw std::runtime_error{"cannot read " + std::string{name}};
    }
    return rows;
}

// A recorded screen: the primers of the Active Primer List ("list",
// comma-separated) and the fields A to D typed on the Multi Primer form
// ("selection": a list number or '_' per field)
struct RecordedSet {
    std::string selection;
    std::string list;

    auto operator<=>(RecordedSet const&) const = default;
};

auto analyse(RecordedSet const& set) -> LegacyMultiPrimer {
    auto listed = std::vector<std::vector<Nucleotide>>{};
    auto stream = std::istringstream{set.list};
    for (auto primer = std::string{}; std::getline(stream, primer, ',');) {
        listed.push_back(nucleotides(primer));
    }
    auto fields = MultiPrimerFields{};
    auto field = std::size_t{0};
    for (auto const typed : set.selection) {
        if (typed != '_') {
            fields.at(field) = listed.at(static_cast<std::size_t>(typed - '1'));
        }
        ++field;
    }
    return legacy_multi_primer(fields, LegacyCriteria{});
}

auto recorded_set(std::vector<std::string> const& row) -> RecordedSet {
    return {.selection = row.at(0), .list = row.at(1)};
}

auto index(std::string const& field) -> std::size_t {
    return static_cast<std::size_t>(field.at(0) - 'A');
}

auto flagged(auto const value, bool const flag) -> std::string {
    return std::to_string(value) + (flag ? "!" : "");
}

// the Balance and Cautions texts of a column, from Length to Hairpins
auto balance_texts(LegacyMultiPrimerColumn const& column) -> std::vector<std::string> {
    return {std::to_string(column.length),
            flagged(column.actual.gc_percent, column.cautions.gc),
            flagged(column.actual.tm_celsius, column.cautions.tm),
            flagged(column.actual.gc_run_3prime, column.cautions.gc_run_3prime),
            flagged(column.actual.longest_run, column.cautions.runs),
            column.cautions.hairpin ? "!" : "OK"};
}

// columns: selection, list, field, position, primer, then length to
// hairpins
void check_balance(std::vector<std::string> const& row) {
    CAPTURE(row);
    constexpr auto first_value = 5;
    auto const column = analyse(recorded_set(row)).columns.at(index(row.at(2)));
    if (not column.has_value()) {
        FAIL("no column for field " << row.at(2));
        return;
    }
    CHECK(balance_texts(*column) == std::vector<std::string>(std::next(row.begin(), first_value), row.end()));
}

// columns: selection, list, row, column, dimer_3prime, dimer_any
void check_dimer(std::vector<std::string> const& row) {
    CAPTURE(row);
    auto const cell = analyse(recorded_set(row)).dimers.at(index(row.at(2))).at(index(row.at(3)));
    if (not cell.has_value()) {
        FAIL("no cell " << row.at(2) << row.at(3));
        return;
    }
    CHECK(flagged(cell->dimer_3prime_matches, cell->flags.three_prime) == row.at(4));
    CHECK(flagged(cell->dimer_adjacent_bases, cell->flags.adjacent) == row.at(5));
}

auto computed_cells(LegacyMultiPrimer const& analysis) -> std::size_t {
    return std::accumulate(analysis.dimers.begin(), analysis.dimers.end(), std::size_t{0},
                           [](std::size_t const sum, auto const& line) {
                               return sum + static_cast<std::size_t>(std::ranges::count_if(
                                                line, [](auto const& cell) { return cell.has_value(); }));
                           });
}

}  // namespace

// Reference values: PRIMER.exe, Analyze -> Multi Primer
// (tests/legacy/multi/multi_primer_balance.tsv and
// multi_primer_dimers.tsv, recorded 2026-10-19)

TEST_CASE("Multi Primer: balance and cautions of every recorded primer", "[analysis][legacy]") {
    for (auto const& row : read_table("multi_primer_balance.tsv")) {
        check_balance(row);
    }
}

TEST_CASE("Multi Primer: every recorded dimer cell, and no other", "[analysis][legacy]") {
    auto recorded_cells = std::map<RecordedSet, std::size_t>{};
    for (auto const& row : read_table("multi_primer_dimers.tsv")) {
        check_dimer(row);
        ++recorded_cells[recorded_set(row)];
    }
    for (auto const& [set, count] : recorded_cells) {
        CAPTURE(set.selection, set.list);
        CHECK(computed_cells(analyse(set)) == count);
    }
}

TEST_CASE("Multi Primer: the diagonal holds self-dimers, the upper triangle pairs", "[analysis][legacy]") {
    // PRIMER.exe (multi_primer_dimers.tsv, selection 1234): A
    // CTCAAAGGCTCAGACATCAG, B GATCCATGCAATGCTAGCTA: A:B 2 and 3, B:B 6!
    auto const analysis = analyse(
        {.selection = "1234",
         .list = "CTCAAAGGCTCAGACATCAG,GATCCATGCAATGCTAGCTA,ACGTTGCAGGTCAAGCTTGC,TTGGCCAAGTCGATCGGATC"});
    auto const& self_b = analysis.dimers.at(1).at(1);
    CHECK((self_b.has_value() and self_b->dimer_3prime_matches == 6 and self_b->flags.three_prime));
    auto const& pair = analysis.dimers.at(0).at(1);
    CHECK((pair.has_value() and pair->dimer_3prime_matches == 2 and pair->dimer_adjacent_bases == 3));
    CHECK(not analysis.dimers.at(1).at(0).has_value());  // shown as dots
}

TEST_CASE("Multi Primer: an empty field leaves its column and cells empty", "[analysis][legacy]") {
    // PRIMER.exe (multi_primer_balance.tsv): fields A = 1, B empty, C = 2
    auto const analysis = analyse({.selection = "1_2", .list = "CTCGATCCGTTCCTAAT,CGGGGCTCATTGCACTC"});
    CHECK(analysis.columns.at(0).has_value());
    CHECK(not analysis.columns.at(1).has_value());
    CHECK(analysis.columns.at(2).has_value());
    CHECK(not analysis.dimers.at(0).at(1).has_value());
    CHECK(analysis.dimers.at(0).at(2).has_value());
}
