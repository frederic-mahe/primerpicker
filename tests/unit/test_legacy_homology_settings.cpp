#include "core/analysis/legacy_homology.hpp"
#include "core/io/genbank.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

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
#include <vector>

namespace {

using primerpicker::analysis::HomologyKind;
using primerpicker::analysis::HomologyOrder;
using primerpicker::analysis::HomologySettings;
using primerpicker::analysis::HomologySite;
using primerpicker::analysis::HomologyTmFilter;
using primerpicker::analysis::legacy_homology_sites;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::Strand;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
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
    auto const parsed =
        primerpicker::io::parse_genbank(read_file(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "homology" / name));
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

auto settings(std::string const& exclude, std::string const& sort) -> HomologySettings {
    return {.tm_filter = exclude == "Y" ? HomologyTmFilter::exclude_at_most_zero : HomologyTmFilter::keep_all,
            .order = sort == "Y" ? HomologyOrder::by_tm : HomologyOrder::by_position};
}

// columns: molecule, primer, position, strand, kind, exclude_tm_at_most_zero,
// sort_by_tm, rank, homologous_bases, site_position, site_strand, tm_celsius
auto recorded_screens() -> std::map<std::vector<std::string>, std::vector<std::string>> {
    constexpr auto column_count = 12U;
    auto screens = std::map<std::vector<std::string>, std::vector<std::string>>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "homology" / "homology_settings.tsv"};
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
        // the rows of a screen are in screen order (rank 1, 2, ...)
        screens[{fields[0], fields[1], fields[4], fields[5], fields[6]}].push_back(
            std::format("{} {} {} {}", fields[8], fields[9], fields[10], fields[11]));
    }
    return screens;
}

}  // namespace

// Reference values: PRIMER.exe, Uniqueness and False Priming screens
// with the System Defaults "exclude Tm <= 0" and "sort homologies by Tm"
// set to Y or N (tests/legacy/homology/homology_settings.tsv, recorded
// 2026-10-19)

TEST_CASE("homologous sites with every recorded setting", "[analysis][legacy]") {
    auto const screens = recorded_screens();
    REQUIRE(screens.size() == 17);
    for (auto const& [key, expected] : screens) {
        INFO("molecule " << key[0] << ", primer " << key[1] << ", " << key[2] << ", exclude " << key[3]
                         << ", sort " << key[4]);
        auto const kind = key[2] == "uniqueness" ? HomologyKind::uniqueness : HomologyKind::false_priming;
        auto const sites = legacy_homology_sites(nucleotides(key[1]), molecule(key[0]), kind, settings(key[3], key[4]));
        auto computed = std::vector<std::string>{};
        std::ranges::transform(sites, std::back_inserter(computed),
                               [&key](HomologySite const& site) { return describe(site, key[1]); });
        CHECK(computed == expected);
    }
}

TEST_CASE("without sorting, false priming sites follow the molecule", "[analysis][legacy]") {
    // PRIMER.exe (homology_settings.tsv): THOM.GB, primer at 51, both
    // settings N: positions 149, 207 (C), 249, 307 (C), 350, 408 (C), 450,
    // 505 (C), whatever the strand
    auto const text = std::string{"GGGGTCGAGCAAGACCGTAT"};
    auto const sites = legacy_homology_sites(
        nucleotides(text), molecule("THOM.GB"), HomologyKind::false_priming,
        {.tm_filter = HomologyTmFilter::keep_all, .order = HomologyOrder::by_position});
    auto positions = std::vector<std::size_t>{};
    std::ranges::transform(sites, std::back_inserter(positions), &HomologySite::position);
    CHECK(positions == std::vector<std::size_t>{149, 207, 249, 307, 350, 408, 450, 505});
}

TEST_CASE("without sorting, uniqueness sites follow each strand", "[analysis][legacy]") {
    // PRIMER.exe (homology_settings.tsv): THOM.GB, primer at 51, sort N:
    // 249 (top), then 207, 307 and 408 (complement)
    auto const text = std::string{"GGGGTCGAGCAAGACCGTAT"};
    auto const sites = legacy_homology_sites(nucleotides(text), molecule("THOM.GB"), HomologyKind::uniqueness,
                                             {.tm_filter = HomologyTmFilter::exclude_at_most_zero,
                                              .order = HomologyOrder::by_position});
    auto positions = std::vector<std::size_t>{};
    std::ranges::transform(sites, std::back_inserter(positions), &HomologySite::position);
    CHECK(positions == std::vector<std::size_t>{249, 207, 307, 408});
}

TEST_CASE("sites with a Tm <= 0 are listed with Tm 0 when kept", "[analysis][legacy]") {
    // PRIMER.exe (homology_settings.tsv): TREP2.GB, primer at 51, exclude
    // N: CATCAG at 294, Tm 0, after the sites of Tm 18 and 10
    auto const text = std::string{"CTCAAAGGCTCAGACATCAG"};
    auto const sites = legacy_homology_sites(nucleotides(text), molecule("TREP2.GB"), HomologyKind::false_priming,
                                             {.tm_filter = HomologyTmFilter::keep_all, .order = HomologyOrder::by_tm});
    REQUIRE(sites.size() == 3);
    CHECK(describe(sites[2], text) == "CATCAG 294 top 0");
}

TEST_CASE("false priming sites of equal Tm keep the molecule order", "[analysis][legacy]") {
    // PRIMER.exe (homology_settings.tsv): THOM3.GB, primer at 51, default
    // settings: CATGCAGAT at 208 (complement) before CATGCAGAT at 400
    // (top), both Tm 10
    auto const text = std::string{"AGCAATTTGATCATGCAGAT"};
    auto const sites = legacy_homology_sites(nucleotides(text), molecule("THOM3.GB"), HomologyKind::false_priming);
    REQUIRE(sites.size() == 2);
    CHECK(describe(sites[0], text) == "CATGCAGAT 208 complement 10");
    CHECK(describe(sites[1], text) == "CATGCAGAT 400 top 10");
}
