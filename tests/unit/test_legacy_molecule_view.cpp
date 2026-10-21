#include "core/io/legacy_molecule_view.hpp"
#include "core/io/molecule_file.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using primerpicker::io::legacy_annotation_lines;
using primerpicker::io::legacy_feature_lines;
using primerpicker::io::legacy_sequence_lines;
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

auto test600_text() -> std::string {
    return read_file(legacy_dir / "pairs" / "TEST600.GB");
}

auto test600() -> std::vector<Nucleotide> {
    auto const parsed = primerpicker::io::parse_molecule_file({.text = test600_text(), .name = "TEST600.GB"});
    auto sequence = parsed.molecule.has_value() ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                                                : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot parse TEST600.GB"};
    }
    return std::move(*sequence);
}

// The lines of the last screen of a recording inside the view's box:
// from the box's second line to the line before its bottom, the text
// from column 2 (after the border and a space) to column 77 (before
// the scroll bar), trailing spaces removed, blank lines left out
auto box_lines(std::filesystem::path const& recording) -> std::vector<std::string> {
    auto const text = read_file(recording);
    auto const last = text.rfind("==== ");
    auto stream = std::istringstream{text.substr(last)};
    auto lines = std::vector<std::string>{};
    auto inside = false;
    for (auto line = std::string{}; std::getline(stream, line);) {
        if (line.starts_with("┌")) {
            inside = true;
            continue;
        }
        // the bottom border, its first cell sometimes under the cursor ('?')
        if (line.starts_with("└") or line.starts_with("?─")) {
            break;
        }
        constexpr auto first_byte = std::size_t{4};  // "│" in UTF-8 (3 bytes), then a space
        constexpr auto content_columns = std::size_t{76};
        if (not inside or line.size() <= first_byte) {
            continue;
        }
        auto content = line.substr(first_byte, content_columns);
        content.erase(content.find_last_not_of(' ') + 1);
        if (not content.empty()) {
            lines.push_back(content);
        }
    }
    return lines;
}

}  // namespace

// Reference: PRIMER.exe, Molecule-Ops -> View/Edit -> Sequence on
// pairs/TEST600.GB (tests/legacy/view_edit/seq_edit.screens.txt, recorded
// 2026-10-21): "      1  GGTCATTTTA AATGGACAAC ..." from column 2
TEST_CASE("View/Edit: the Sequence view's lines", "[io][view][legacy]") {
    auto const recorded = box_lines(legacy_dir / "view_edit" / "seq_edit.screens.txt");
    REQUIRE(recorded.size() == 10);
    CHECK(legacy_sequence_lines(test600()) == recorded);
}

// Reference: the same, View/Edit -> Annotations
// (view_edit/annotations.screens.txt): the file's lines up to ORIGIN
TEST_CASE("View/Edit: the Annotations view shows the file's header lines", "[io][view][legacy]") {
    auto const recorded = box_lines(legacy_dir / "view_edit" / "annotations.screens.txt");
    REQUIRE(recorded.size() == 4);
    CHECK(legacy_annotation_lines(test600_text()) == recorded);
}

// Reference: the same, View/Edit -> Features
// (view_edit/features.screens.txt): "-- no features --"
TEST_CASE("View/Edit: a molecule without features", "[io][view][legacy]") {
    CHECK(legacy_feature_lines(test600_text()) == box_lines(legacy_dir / "view_edit" / "features.screens.txt"));
}

// Not recorded (no molecule with features was probed): the lines of
// the FEATURES block, as written in the file (GenBank feature table)
TEST_CASE("View/Edit: the lines of a features table", "[io][view]") {
    constexpr auto text = std::string_view{
        "LOCUS       TFEAT          30 bp    DNA     linear\r\n"
        "FEATURES             Location/Qualifiers\r\n"
        "     gene            3..20\r\n"
        "                     /gene=\"abc\"\r\n"
        "ORIGIN\r\n"
        "        1 acgtacgtac gtacgtacgt acgtacgtac\r\n"
        "//\r\n"};
    CHECK(legacy_feature_lines(text) ==
          std::vector<std::string>{"     gene            3..20", "                     /gene=\"abc\""});
}
