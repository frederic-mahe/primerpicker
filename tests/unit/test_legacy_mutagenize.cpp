#include "core/analysis/legacy_mutagenize.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/io/genbank.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <charconv>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <iterator>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using primerpicker::analysis::legacy_mutagenize_view;
using primerpicker::sequence::Nucleotide;

auto read_file(std::filesystem::path const& path) -> std::string {
    auto text = std::string(std::filesystem::file_size(path), '\0');
    auto file = std::ifstream{path, std::ios::binary};
    file.read(text.data(), std::ssize(text));
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    return text;
}

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

auto test600() -> std::vector<Nucleotide> {
    auto const parsed =
        primerpicker::io::parse_genbank(read_file(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "pairs" / "TEST600.GB"));
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"cannot read TEST600.GB"};
    }
    return nucleotides(parsed.molecule->sequence);
}

auto number(std::string const& text) -> std::size_t {
    auto value = std::size_t{0};
    auto const* const last = std::next(text.data(), std::ssize(text));
    auto const [end, error] = std::from_chars(text.data(), last, value);
    if (error != std::errc{} or end != last) {
        throw std::runtime_error{"not a number: " + text};
    }
    return value;
}

// mutagenize/mutagenize.tsv: columns molecule, position, keys, line,
// text; one screen per (position, keys)
using Screen = std::map<std::string, std::string>;

auto recorded_screens() -> std::map<std::pair<std::size_t, std::string>, Screen> {
    constexpr auto column_count = 5U;
    auto screens = std::map<std::pair<std::size_t, std::string>, Screen>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "mutagenize" / "mutagenize.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        // a line may end with an empty text
        if (fields.size() == column_count - 1) {
            fields.emplace_back();
        }
        if (fields.size() != column_count or fields[0] != "TEST600.GB") {
            throw std::runtime_error{"unexpected line: " + line};
        }
        screens[{number(fields[1]), fields[2]}][fields[3]] = fields[4];
    }
    return screens;
}

// the letters of the primer line (the mutated primer)
auto primer_letters(std::string const& line) -> std::string {
    auto letters = std::string{};
    for (auto const letter : line) {
        if (letter != ' ') {
            letters += letter;
        }
    }
    return letters;
}

void check_screen(std::vector<Nucleotide> const& molecule, std::size_t const position, Screen const& screen) {
    CAPTURE(position);
    auto const mutated = nucleotides(primer_letters(screen.at("primer")));
    auto const view = legacy_mutagenize_view(molecule, position, mutated);
    auto const computed = std::vector<std::string>{view.original_frames.at(0), view.original_frames.at(1),
                                                   view.original_frames.at(2), view.molecule,
                                                   view.primer, view.mutated_frames.at(0),
                                                   view.mutated_frames.at(1), view.mutated_frames.at(2)};
    auto expected = std::vector<std::string>{};
    for (auto const* const name : {"O", "R", "I", "molecule", "primer", "M", "U", "T"}) {
        expected.push_back(screen.at(name));
    }
    CHECK(computed == expected);
}

// "%GC: 50    Tm°C: 66    Hairpin: OK  Dimers: OK  Runs: !"
auto summary(std::vector<Nucleotide> const& primer) -> std::string {
    auto const actual = primerpicker::analysis::legacy_actual_values(primer);
    auto const verdict = primerpicker::analysis::evaluate(actual, primerpicker::analysis::LegacyCriteria{});
    auto const flag = [](bool const meets) { return meets ? std::string{"OK"} : std::string{"!"}; };
    return std::format("%GC: {:<6}Tm°C: {:<6}Hairpin: {:<4}Dimers: {:<4}Runs: {}",
                       std::to_string(actual.gc_percent) + (verdict.gc ? "" : "!"),
                       std::to_string(actual.tm_celsius) + (verdict.tm ? "" : "!"), flag(verdict.no_hairpin),
                       flag(verdict.no_3prime_dimer and verdict.no_dimer),
                       flag(verdict.no_runs and verdict.no_3prime_gc_run));
}

}  // namespace

// Reference values: PRIMER.exe, Modify -> Mutagenize on TEST600.GB
// (tests/legacy/mutagenize/mutagenize.tsv, recorded 2026-10-19)

TEST_CASE("Mutagenize: every recorded screen", "[analysis][legacy]") {
    auto const molecule = test600();
    auto const screens = recorded_screens();
    REQUIRE(screens.size() == 6);
    for (auto const& [key, screen] : screens) {
        check_screen(molecule, key.first, screen);
    }
}

TEST_CASE("Mutagenize: the values of the original and mutated primers", "[analysis][legacy]") {
    auto const molecule = test600();
    for (auto const& [key, screen] : recorded_screens()) {
        CAPTURE(key.first, key.second);
        auto const primer_length = primer_letters(screen.at("primer")).size();
        auto const first = std::next(molecule.begin(), static_cast<std::ptrdiff_t>(key.first - 1));
        auto const original = std::vector<Nucleotide>(first, std::next(first, static_cast<std::ptrdiff_t>(primer_length)));
        CHECK(summary(original) == screen.at("original"));
        CHECK(summary(nucleotides(primer_letters(screen.at("primer")))) == screen.at("mutated"));
    }
}

TEST_CASE("Mutagenize: amino acids under the third base of their codon", "[analysis][legacy]") {
    // PRIMER.exe (mutagenize.tsv): at position 575 the window ends with
    // the molecule (600 bases); codon 598-600 GAC is D, under base 600
    auto const molecule = test600();
    auto const primer = std::vector<Nucleotide>(std::next(molecule.begin(), 574), std::next(molecule.begin(), 594));
    auto const view = legacy_mutagenize_view(molecule, 575, primer);
    CHECK(view.original_frames.at(2).ends_with("S  D"));
    CHECK(view.molecule.ends_with("GCGAC"));
}
