#include "core/analysis/legacy_homology.hpp"
#include "core/io/genbank.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
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
#include <utility>
#include <vector>

namespace {

using primerpicker::analysis::HomologyKind;
using primerpicker::analysis::HomologySite;
using primerpicker::analysis::legacy_homology_sites;
using primerpicker::analysis::Strand;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::to_nucleotides;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

auto read_file(std::filesystem::path const& path) -> std::string {
    auto text = std::string(std::filesystem::file_size(path), '\0');
    auto file = std::ifstream{path, std::ios::binary};
    file.read(text.data(), std::ssize(text));
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    return text;
}

auto molecule(std::string const& name) -> std::vector<Nucleotide> {
    auto const directory = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / (name.starts_with("TREP") ? "homology" : "pairs");
    auto const parsed = primerpicker::io::parse_genbank(read_file(directory / name));
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"cannot parse " + name};
    }
    return nucleotides(parsed.molecule->sequence);
}

// one site as recorded: bases, position, strand, Tm
auto describe(HomologySite const& site, std::string_view const primer) -> std::string {
    return std::format("{} {} {} {}", primer.substr(site.primer_start, site.length), site.position,
                       site.strand == Strand::complement ? "complement" : "top", site.tm_celsius);
}

}  // namespace

// Reference values: PRIMER.exe, Uniqueness and False Priming screens
// (tests/legacy/homology/homology.tsv, recorded 2026-10-18)

TEST_CASE("a unique primer has no homologous site", "[analysis][legacy]") {
    auto const template_sequence = molecule("TEST600.GB");
    auto const primer = nucleotides("CTCAAAGGCTCAGACATCAG");
    CHECK(legacy_homology_sites(primer, template_sequence, HomologyKind::uniqueness).empty());  // PRIMER.exe
}

TEST_CASE("false priming sites include the 3' end", "[analysis][legacy]") {
    auto const template_sequence = molecule("TREP2.GB");
    auto const text = std::string{"CTCAAAGGCTCAGACATCAG"};
    auto const sites = legacy_homology_sites(nucleotides(text), template_sequence, HomologyKind::false_priming);
    REQUIRE(sites.size() == 2);
    CHECK(describe(sites[0], text) == "TCAGACATCAG 230 complement 18");  // PRIMER.exe
    CHECK(describe(sites[1], text) == "CAGACATCAG 150 top 10");         // PRIMER.exe
}

TEST_CASE("homologous sites match the recorded screens", "[analysis][legacy]") {
    // key: molecule, primer, kind; value: the recorded sites, in order
    auto recorded = std::map<std::vector<std::string>, std::vector<std::string>>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "homology" / "homology.tsv"};
    constexpr auto column_count = 9U;
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
        // cppcheck-suppress variableScope ; operator[] creates the key, even without a site
        auto& sites = recorded[{fields[0], fields[1], fields[4]}];
        if (fields[5] != "-") {
            sites.push_back(std::format("{} {} {} {}", fields[5], fields[6], fields[7], fields[8]));
        }
    }
    REQUIRE(recorded.size() == 7);
    for (auto const& [key, expected] : recorded) {
        INFO("molecule " << key[0] << ", primer " << key[1] << ", " << key[2]);
        auto const kind = key[2] == "uniqueness" ? HomologyKind::uniqueness : HomologyKind::false_priming;
        auto const sites = legacy_homology_sites(nucleotides(key[1]), molecule(key[0]), kind);
        auto computed = std::vector<std::string>{};
        std::ranges::transform(sites, std::back_inserter(computed),
                               [&key](HomologySite const& site) { return describe(site, key[1]); });
        CHECK(computed == expected);
    }
}
