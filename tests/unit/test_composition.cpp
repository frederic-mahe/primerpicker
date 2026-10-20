#include "core/sequence/composition.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using primerpicker::sequence::Composition;
using primerpicker::sequence::composition;
using primerpicker::sequence::Nucleotide;

auto const reports_dir = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "reports";

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

auto read_file(std::filesystem::path const& path) -> std::string {
    auto file = std::ifstream{path, std::ios::binary};
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    auto text = std::ostringstream{};
    text << file.rdbuf();
    return std::move(text).str();
}

// "Composition:    A: 6  C: 4  G: 5  T: 2", the line of a One Primer
// report, without its indentation and line end
auto composition_line(std::string_view const report) -> std::string {
    constexpr auto label = std::string_view{"Composition:"};
    auto const start = report.find(label);
    if (start == std::string_view::npos) {
        throw std::runtime_error{"no Composition line"};
    }
    auto const end = report.find('\r', start);
    return std::string{report.substr(start, end - start)};
}

auto shown(Composition const& counts) -> std::string {
    return std::format("Composition:    A:{:2}  C:{:2}  G:{:2}  T:{:2}", counts.adenine, counts.cytosine,
                       counts.guanine, counts.thymine);
}

}  // namespace

// Reference: PRIMER.exe, One Primer report of the M13 primer
// (tests/legacy/reports/INSPECT.PRN, recorded 2026-10-19)
TEST_CASE("composition: counts of each nucleotide", "[composition]") {
    auto const primer = nucleotides("GTAAAACGACGGCCAGT");
    auto const counts = composition(primer);
    CHECK(counts.adenine == 6);
    CHECK(counts.cytosine == 4);
    CHECK(counts.guanine == 5);
    CHECK(counts.thymine == 2);
}

TEST_CASE("composition of an empty sequence", "[composition]") {
    STATIC_REQUIRE(composition({}) == Composition{});
}

// Reference: the Composition line of every recorded One Primer report
// (tests/legacy/reports/inspect_reports.tsv, recorded 2026-10-20); it
// counts the 5' extension too
TEST_CASE("composition: every recorded One Primer report", "[composition][legacy]") {
    auto table = std::ifstream{reports_dir / "inspect_reports.tsv"};
    REQUIRE(table);
    auto checked = std::size_t{0};
    for (auto row = std::string{}; std::getline(table, row);) {
        if (row.starts_with('#')) {
            continue;
        }
        auto fields = std::istringstream{row};
        auto file = std::string{};
        auto primer = std::string{};
        auto extension = std::string{};
        fields >> file >> primer >> extension;
        auto const sequence = nucleotides((extension == "-" ? std::string{} : extension) + primer);
        INFO(file);
        CHECK(shown(composition(sequence)) == composition_line(read_file(reports_dir / file)));
        ++checked;
    }
    CHECK(checked == 47);
}
