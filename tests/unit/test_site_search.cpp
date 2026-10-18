#include "core/enzymes/site_search.hpp"
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
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using primerpicker::enzymes::find_sites;
using primerpicker::enzymes::Mismatches;
using primerpicker::enzymes::parse_site;
using primerpicker::enzymes::Site;
using primerpicker::enzymes::SitePattern;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::Strand;

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
    return std::move(*converted);
}

auto tenz() -> std::vector<Nucleotide> {
    auto const path = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "enzymes" / "TENZ.GB";
    auto const parsed = primerpicker::io::parse_genbank(read_file(path));
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"cannot read TENZ.GB"};
    }
    return nucleotides(parsed.molecule->sequence);
}

auto pattern(std::string_view const site) -> SitePattern {
    auto parsed = parse_site(site);
    if (not parsed.has_value()) {
        throw std::runtime_error{"not a recognition sequence: " + std::string{site}};
    }
    return std::move(*parsed);
}

auto letter(Nucleotide const nucleotide) -> char {
    constexpr auto letters = std::string_view{"ACGT"};
    return letters.at(static_cast<std::size_t>(nucleotide));
}

// a site as the recorded tables write it: position, strand, and for
// almost sites a dot per matching base, the molecule's base otherwise
auto describe(Site const& site, std::vector<Nucleotide> const& molecule, std::size_t const length,
              Mismatches const mismatches) -> std::string {
    auto const strand = std::string_view{site.strand == Strand::top ? "top" : "complement"};
    if (mismatches == Mismatches::none) {
        return std::format("{} {} -", site.position, strand);
    }
    auto match = std::string(length, '.');
    if (site.mismatch_offset.has_value()) {
        match.at(*site.mismatch_offset) = letter(molecule.at(site.position - 1 + *site.mismatch_offset));
    }
    return std::format("{} {} {}", site.position, strand, match);
}

auto split(std::string const& line) -> std::vector<std::string> {
    auto fields = std::vector<std::string>{};
    auto stream = std::istringstream{line};
    for (auto field = std::string{}; std::getline(stream, field, '\t');) {
        fields.push_back(field);
    }
    return fields;
}

struct RecordedSearch {
    std::string kind;
    std::string enzyme;
    std::string site;

    friend auto operator<=>(RecordedSearch const&, RecordedSearch const&) = default;
};

// tests/legacy/enzymes/find_enzymes.tsv, by search: sites by rank
auto recorded_searches() -> std::map<RecordedSearch, std::vector<std::string>> {
    auto searches = std::map<RecordedSearch, std::vector<std::string>>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "enzymes" / "find_enzymes.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto const fields = split(line);
        searches[{.kind = fields.at(1), .enzyme = fields.at(2), .site = fields.at(3)}].push_back(
            std::format("{} {} {}", fields.at(5), fields.at(6), fields.at(7)));
    }
    return searches;
}

}  // namespace

// Reference values: PRIMER.exe, True Sites and Almost Sites screens on
// TENZ.GB (tests/legacy/enzymes/find_enzymes.tsv, recorded 2026-10-18)

TEST_CASE("Site search finds the recorded true and almost sites", "[enzymes][legacy]") {
    auto const molecule = tenz();
    auto const searches = recorded_searches();
    REQUIRE(searches.size() == 6);
    for (auto const& [search, expected] : searches) {
        auto const mismatches = search.kind == "almost" ? Mismatches::at_most_one : Mismatches::none;
        auto const sites = find_sites(molecule, pattern(search.site), mismatches);
        auto actual = std::vector<std::string>{};
        std::ranges::transform(sites, std::back_inserter(actual), [&](Site const& site) {
            return describe(site, molecule, search.site.size(), mismatches);
        });
        INFO(search.kind << " " << search.enzyme);
        CHECK(actual == expected);
    }
}

// The cases below follow from the IUPAC codes (NC-IUB 1985) and from
// the rules of the recorded screens (tests/legacy/README.md, "True and
// almost sites")

TEST_CASE("Site search: IUPAC letters in the recognition sequence", "[enzymes]") {
    // AccI GTMKAC (M: A or C; K: G or T), a palindrome: GTATAC and GTCGAC
    // match, GTGTAC does not
    auto const molecule = nucleotides("GTATACTTGTCGACTTGTGTAC");
    auto const sites = find_sites(molecule, pattern("GTMKAC"), Mismatches::none);
    REQUIRE(sites.size() == 2);
    CHECK(sites[0] == Site{.position = 1, .strand = Strand::top, .mismatch_offset = std::nullopt});
    CHECK(sites[1] == Site{.position = 9, .strand = Strand::top, .mismatch_offset = std::nullopt});
}

TEST_CASE("Site search: N matches any base", "[enzymes]") {
    // BglI GCCNNNNNGGC (palindrome)
    auto const molecule = nucleotides("TTGCCATGCAGGCTT");
    auto const sites = find_sites(molecule, pattern("GCCNNNNNGGC"), Mismatches::none);
    REQUIRE(sites.size() == 1);
    CHECK(sites[0].position == 3);
}

TEST_CASE("Site search: non-palindromic sites on the complement strand", "[enzymes]") {
    // AlwI GGATC; its reverse complement GATCC at 6 is a complement site
    auto const molecule = nucleotides("GGATCGATCC");
    auto const sites = find_sites(molecule, pattern("GGATC"), Mismatches::none);
    REQUIRE(sites.size() == 2);
    CHECK(sites[0] == Site{.position = 1, .strand = Strand::top, .mismatch_offset = std::nullopt});
    CHECK(sites[1] == Site{.position = 6, .strand = Strand::complement, .mismatch_offset = std::nullopt});
}

TEST_CASE("Site search: one mismatch, and its offset", "[enzymes]") {
    // ApaI GGGCCC: GGGCCA differs at offset 5
    auto const molecule = nucleotides("TTGGGCCATT");
    CHECK(find_sites(molecule, pattern("GGGCCC"), Mismatches::none).empty());
    auto const sites = find_sites(molecule, pattern("GGGCCC"), Mismatches::at_most_one);
    REQUIRE(sites.size() == 1);
    CHECK(sites[0] == Site{.position = 3, .strand = Strand::top, .mismatch_offset = 5});
}

TEST_CASE("Site search: a molecule shorter than the site has no site", "[enzymes]") {
    CHECK(find_sites(nucleotides("GAATT"), pattern("GAATTC"), Mismatches::at_most_one).empty());
    CHECK(find_sites(std::vector<Nucleotide>{}, pattern("GAATTC"), Mismatches::none).empty());
}

TEST_CASE("Recognition sequences: IUPAC letters only", "[enzymes]") {
    CHECK(parse_site("GAATTC").has_value());
    CHECK(parse_site("gcnngc").has_value());
    CHECK(not parse_site("").has_value());
    CHECK(not parse_site("GAAUUC").has_value());
    CHECK(not parse_site("GAA TTC").has_value());
}
