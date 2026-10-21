#include "core/analysis/legacy_verdict.hpp"
#include "core/io/molecule_file.hpp"
#include "core/search/legacy_picking.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using primerpicker::analysis::LegacyCriteria;
using primerpicker::search::legacy_mutagenize_line;
using primerpicker::sequence::Nucleotide;

auto const legacy_dir = std::filesystem::path{PRIMERPICKER_LEGACY_DIR};

auto read_file(std::filesystem::path const& path) -> std::string {
    auto file = std::ifstream{path, std::ios::binary};
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    auto text = std::ostringstream{};
    text << file.rdbuf();
    return std::move(text).str();
}

auto molecule(std::filesystem::path const& relative_path) -> std::vector<Nucleotide> {
    auto const text = read_file(legacy_dir / relative_path);
    auto const parsed =
        primerpicker::io::parse_molecule_file({.text = text, .name = relative_path.filename().string()});
    auto sequence = parsed.molecule.has_value() ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                                                : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot parse " + relative_path.string()};
    }
    return std::move(*sequence);
}

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

// the bases of the primer line of the screen
auto primer_letters(std::string const& line) -> std::string {
    auto letters = std::string{};
    for (auto const letter : line) {
        if (letter != ' ') {
            letters += letter;
        }
    }
    return letters;
}

// One recorded screen: its primer position, the primer line, and the
// ORIGINAL and MUTATED lines
struct RecordedLines {
    std::size_t position = 0;
    std::string primer;
    std::string original;
    std::string mutated;
};

// mutagenize/mutagenize.tsv: molecule, position, keys, line, text
auto recorded_lines() -> std::vector<RecordedLines> {
    auto screens = std::map<std::pair<std::string, std::string>, RecordedLines>{};
    auto file = std::ifstream{legacy_dir / "mutagenize" / "mutagenize.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        constexpr auto text_column = std::size_t{4};
        if (fields.size() <= text_column) {
            continue;  // an empty translation row
        }
        auto& screen = screens[{fields[1], fields[2]}];
        screen.position = std::stoul(fields[1]);
        if (fields[3] == "primer") {
            screen.primer = primer_letters(fields[text_column]);
        } else if (fields[3] == "original") {
            screen.original = fields[text_column];
        } else if (fields[3] == "mutated") {
            screen.mutated = fields[text_column];
        }
    }
    auto found = std::vector<RecordedLines>{};
    for (auto const& [key, screen] : screens) {
        found.push_back(screen);
    }
    return found;
}

void check_lines(std::vector<Nucleotide> const& sequence, RecordedLines const& recorded) {
    INFO(recorded.position << " " << recorded.primer);
    auto const mutated = nucleotides(recorded.primer);
    auto const first = std::next(sequence.begin(), static_cast<std::ptrdiff_t>(recorded.position - 1));
    auto const original = std::vector<Nucleotide>(first, std::next(first, std::ssize(mutated)));
    CHECK(legacy_mutagenize_line(original, LegacyCriteria{}) == recorded.original);
    CHECK(legacy_mutagenize_line(mutated, LegacyCriteria{}) == recorded.mutated);
}

}  // namespace

// Reference: PRIMER.exe, Modify -> Mutagenize on TEST600.GB, the
// ORIGINAL and MUTATED lines of every recorded screen
// (tests/legacy/mutagenize/mutagenize.tsv, recorded 2026-10-19)
TEST_CASE("Mutagenize: the ORIGINAL and MUTATED lines of every recorded screen", "[search][mutagenize][legacy]") {
    auto const sequence = molecule("pairs/TEST600.GB");
    auto const screens = recorded_lines();
    REQUIRE(screens.size() == 6);
    for (auto const& recorded : screens) {
        check_lines(sequence, recorded);
    }
}

// Reference: PRIMER.exe, Mutagenize on the primer at 450 of
// search/TRUN.GB, which ends with AAA (tests/legacy/mutagenize_3prime_run/
// TRUN_450.screens.txt, recorded 2026-10-21): the run at the 3' end
// counts, as on the picking screen
TEST_CASE("Mutagenize: the run at the 3' end counts", "[search][mutagenize][legacy]") {
    CHECK(legacy_mutagenize_line(nucleotides("GCAGTCAGGCTCATCGCAAA"), LegacyCriteria{}) ==
          "%GC: 55    Tm°C: 75    Hairpin: OK  Dimers: OK  Runs: !");
}
