#include "core/analysis/legacy_hairpin.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using Catch::Matchers::WithinAbs;
using primerpicker::analysis::legacy_hairpin;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::to_nucleotides;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

// throws (and fails the test) when there is no hairpin
auto hairpin_of(std::string_view const primer) -> primerpicker::analysis::Hairpin {
    auto const hairpin = legacy_hairpin(nucleotides(primer));
    if (not hairpin.has_value()) {
        throw std::runtime_error{"no hairpin: " + std::string{primer}};
    }
    return *hairpin;
}

// energies are shown with one decimal
constexpr auto tolerance = 0.05;

}  // namespace

// Reference values: PRIMER.exe, Direct Entry -> Inspect -> HairPn
// (tests/legacy/inspect/), built-in defaults

TEST_CASE("hairpin of a designed stem-loop", "[analysis][legacy]") {
    auto const hairpin = hairpin_of("GACGCTTTTGCGTC");
    CHECK_THAT(hairpin.energy.kcal_per_mol, WithinAbs(-3.9, tolerance));  // PRIMER.exe
    CHECK(hairpin.stem_length == 5);  // PRIMER.exe: ":::::"
    CHECK(hairpin.loop_length == 4);  // PRIMER.exe
    CHECK(hairpin.stem_start == 0);
}

TEST_CASE("hairpin energy depends on the loop length", "[analysis][legacy]") {
    auto const energy = [](std::string_view const primer) {
        auto const hairpin = legacy_hairpin(nucleotides(primer));
        return hairpin.has_value() ? hairpin->energy.kcal_per_mol : 0.0;
    };
    CHECK_THAT(energy("GACGCTTTGCGTC"), WithinAbs(-2.4, tolerance));         // PRIMER.exe, loop 3
    CHECK_THAT(energy("GACGCTTTTTGCGTC"), WithinAbs(-5.4, tolerance));       // PRIMER.exe, loop 5
    CHECK_THAT(energy("GACGCTTTTTTTTTTGCGTC"), WithinAbs(-5.5, tolerance));  // PRIMER.exe, loop 10
}

TEST_CASE("the loop is zipped down to four bases in repeats", "[analysis][legacy]") {
    auto const hairpin = hairpin_of("ATATATATATATATATATAT");
    CHECK_THAT(hairpin.energy.kcal_per_mol, WithinAbs(-1.0, tolerance));  // PRIMER.exe
    CHECK(hairpin.loop_length == 4);  // PRIMER.exe
    CHECK(hairpin.stem_length == 8);  // PRIMER.exe
}

TEST_CASE("no hairpin when no stem-loop is stable", "[analysis][legacy]") {
    CHECK_FALSE(legacy_hairpin(nucleotides("GACTTTTGTC")).has_value());  // PRIMER.exe: No Hairpins Found
    CHECK_FALSE(legacy_hairpin(nucleotides("AAAAAAAAAAAAAAAAAAAAAAAA")).has_value());  // PRIMER.exe
}

namespace {

struct Recording {
    std::string primer;
    std::string energy;   // "-" when no hairpin is found
    std::string pairing;  // the marks line of the drawing
    std::string top;      // the 5' part line of the drawing
};

auto read_recordings(std::filesystem::path const& path) -> std::vector<Recording> {
    constexpr auto column_count = 10U;
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
        recordings.push_back({.primer = fields[0], .energy = fields[6], .pairing = fields[8], .top = fields[7]});
    }
    return recordings;
}

auto shown_energy(Recording const& recording) -> std::optional<double> {
    if (recording.energy == "-") {
        return std::nullopt;
    }
    auto stream = std::istringstream{recording.energy};
    auto value = 0.0;
    stream >> value;
    if (stream.fail() or not stream.eof()) {
        throw std::runtime_error{"unexpected energy: " + recording.energy};
    }
    return value;
}

auto computed_energy(std::string const& primer) -> std::optional<double> {
    auto const hairpin = legacy_hairpin(nucleotides(primer));
    if (not hairpin.has_value()) {
        return std::nullopt;
    }
    return hairpin->energy.kcal_per_mol;
}

auto same(std::optional<double> const lhs, std::optional<double> const rhs) -> bool {
    if (lhs.has_value() != rhs.has_value()) {
        return false;
    }
    return not lhs.has_value() or (*lhs - *rhs < tolerance and *rhs - *lhs < tolerance);
}

}  // namespace

TEST_CASE("legacy hairpins match the recorded primers", "[analysis][legacy]") {
    auto const directory = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "inspect";
    auto recordings = read_recordings(directory / "inspect.tsv");
    auto const designed = read_recordings(directory / "hairpin_design.tsv");
    recordings.insert(recordings.end(), designed.begin(), designed.end());
    REQUIRE(recordings.size() > 140);
    for (auto const& recording : recordings) {
        INFO("primer " << recording.primer << ", shown " << recording.energy);
        CHECK(same(computed_energy(recording.primer), shown_energy(recording)));
    }
}

TEST_CASE("hairpin search works at compile time", "[analysis][legacy][constexpr]") {
    // GACGCTTTTGCGTC: -3.9 kcal (PRIMER.exe, hairpin_design.tsv)
    constexpr auto energy = [] {
        auto const hairpin = legacy_hairpin(to_nucleotides("GACGCTTTTGCGTC").value());
        return hairpin.has_value() ? hairpin->energy.kcal_per_mol : 0.0;
    }();
    STATIC_CHECK(energy < -3.85);
    STATIC_CHECK(energy > -3.95);
}
