#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_pair.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/report/legacy_report.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <chrono>
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

using primerpicker::analysis::LegacyCriteria;
using primerpicker::analysis::LegacyMultiPrimerFields;
using primerpicker::analysis::LegacyPrimer;
using primerpicker::analysis::LegacyPrimerPair;
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

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence"};
    }
    return *converted;
}

constexpr auto printed = std::chrono::year_month_day{std::chrono::year{2026}, std::chrono::October, std::chrono::day{22}};

}  // namespace

// Reference: PRIMER.exe, the Primer Pair and Multiple Primers reports of
// M13EXT.PD2 (GAATTC + GTAAAACGACGGCCAGT) and AACAGCTATGACCATG, printed to
// disk on 10/22/26 (tests/legacy/extension/pair_ext_report.PRN,
// multi_ext_report.PRN, recorded 2026-10-22): the sequence with a dash
// after the extension, the values of the whole primer (the Tm without
// the extension), the duplexes of the whole primers
TEST_CASE("Primer Pair report of a primer with a 5' extension", "[report][legacy]") {
    auto const primer_a = nucleotides("GAATTCGTAAAACGACGGCCAGT");
    auto const primer_b = nucleotides("AACAGCTATGACCATG");
    auto const report = primerpicker::report::legacy_pair_report(
        LegacyPrimerPair{.first = {.sequence = primer_a, .extension_length = 6},
                         .second = {.sequence = primer_b, .extension_length = 0}},
        LegacyCriteria{}, printed);
    CHECK(report == read_file(legacy_dir / "extension" / "pair_ext_report.PRN"));
}

TEST_CASE("Multiple Primers report of a primer with a 5' extension", "[report][legacy]") {
    auto const primer_a = nucleotides("GAATTCGTAAAACGACGGCCAGT");
    auto const primer_b = nucleotides("AACAGCTATGACCATG");
    auto const fields =
        LegacyMultiPrimerFields{LegacyPrimer{.sequence = primer_a, .extension_length = 6},
                                LegacyPrimer{.sequence = primer_b, .extension_length = 0}, std::nullopt, std::nullopt};
    auto const report = primerpicker::report::legacy_multi_primer_report(fields, LegacyCriteria{}, printed);
    CHECK(report == read_file(legacy_dir / "extension" / "multi_ext_report.PRN"));
}
