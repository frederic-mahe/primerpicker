#include "core/analysis/legacy_mutagenize.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/io/genbank.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <catch2/catch_test_macros.hpp>

#include <charconv>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using primerpicker::analysis::legacy_mutagenize_strand;
using primerpicker::analysis::legacy_mutagenize_view;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::PrimerLocation;
using primerpicker::sequence::Strand;

constexpr auto primer_length = std::size_t{20};  // the recorded primers

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

// "400C": position 400 of the complement strand
auto complement_position(std::string const& text) -> std::size_t {
    if (not text.ends_with('C')) {
        throw std::runtime_error{"not a complement-strand position: " + text};
    }
    auto value = std::size_t{0};
    auto const* const last = std::next(text.data(), std::ssize(text) - 1);
    auto const [end, error] = std::from_chars(text.data(), last, value);
    if (error != std::errc{} or end != last) {
        throw std::runtime_error{"not a number: " + text};
    }
    return value;
}

// complement/mutagenize_complement.tsv: columns molecule, position,
// keys, line, text; one screen per (position, keys)
using Screen = std::map<std::string, std::string>;

auto recorded_screens() -> std::map<std::pair<std::size_t, std::string>, Screen> {
    constexpr auto column_count = 5U;
    auto screens = std::map<std::pair<std::size_t, std::string>, Screen>{};
    auto file =
        std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "complement" / "mutagenize_complement.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        if (fields.size() != column_count or fields[0] != "TEST600.GB") {
            throw std::runtime_error{"unexpected line: " + line};
        }
        screens[{complement_position(fields[1]), fields[2]}][fields[3]] = fields[4];
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

// Reference values: PRIMER.exe, Modify -> Mutagenize on complement-strand
// primers of TEST600.GB (tests/legacy/complement/mutagenize_complement.tsv,
// recorded 2026-10-22)

TEST_CASE("Mutagenize, complement strand: every recorded screen", "[analysis][legacy]") {
    auto const molecule = test600();
    auto const screens = recorded_screens();
    REQUIRE(screens.size() == 3);
    for (auto const& [key, screen] : screens) {
        CAPTURE(key.first, key.second);
        auto const strand = legacy_mutagenize_strand(
            molecule, {.position = key.first, .length = primer_length, .strand = Strand::complement});
        auto const view =
            legacy_mutagenize_view(strand.molecule, strand.position, nucleotides(primer_letters(screen.at("primer"))));
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
}

TEST_CASE("Mutagenize, complement strand: the original primer is the primer's own strand", "[analysis][legacy]") {
    auto const molecule = test600();
    for (auto const& [key, screen] : recorded_screens()) {
        CAPTURE(key.first, key.second);
        auto const location = PrimerLocation{.position = key.first, .length = primer_length, .strand = Strand::complement};
        auto const strand = legacy_mutagenize_strand(molecule, location);
        auto const first = std::next(strand.molecule.begin(), static_cast<std::ptrdiff_t>(strand.position - 1));
        auto const shown = std::vector<Nucleotide>(first, std::next(first, static_cast<std::ptrdiff_t>(primer_length)));
        CHECK(shown == primerpicker::sequence::primer_at(molecule, location));
        CHECK(summary(shown) == screen.at("original"));
        CHECK(summary(nucleotides(primer_letters(screen.at("primer")))) == screen.at("mutated"));
    }
}

TEST_CASE("Mutagenize, complement strand: the reverse complement of the molecule", "[analysis][legacy]") {
    // PRIMER.exe (complement/mut_c400.screens.txt): for 400(C) of the 600
    // bases of TEST600, the box shows TCATCTGGCG... from 20 bases before
    // the primer's 5' end, which is base 600 - 400 + 1 = 201 of the
    // reverse complement
    auto const molecule = test600();
    auto const strand =
        legacy_mutagenize_strand(molecule, {.position = 400, .length = primer_length, .strand = Strand::complement});
    CHECK(strand.position == 201);
    REQUIRE(strand.molecule.size() == molecule.size());
    CHECK(strand.molecule.front() == primerpicker::sequence::complement(molecule.back()));
    auto const top = legacy_mutagenize_strand(molecule, {.position = 51, .length = primer_length, .strand = Strand::top});
    CHECK(top.position == 51);
    CHECK(top.molecule == molecule);
}
