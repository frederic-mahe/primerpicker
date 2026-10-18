#include "core/enzymes/site_search.hpp"
#include "core/io/genbank.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "core/sequence/topology.hpp"

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
using primerpicker::sequence::Topology;

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

struct RecordedMolecule {
    std::vector<Nucleotide> sequence;
    Topology topology;
};

auto recorded_molecule(std::string const& name) -> RecordedMolecule {
    auto const path = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "enzymes" / name;
    auto const parsed = primerpicker::io::parse_genbank(read_file(path));
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"cannot read " + name};
    }
    return {.sequence = nucleotides(parsed.molecule->sequence), .topology = parsed.molecule->topology};
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
    if (site.mismatch_offset.has_value()) {  // circular molecules: across the origin
        match.at(*site.mismatch_offset) =
            letter(molecule.at((site.position - 1 + *site.mismatch_offset) % molecule.size()));
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
    std::string molecule;
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
        auto const search =
            RecordedSearch{.molecule = fields.at(0), .kind = fields.at(1), .enzyme = fields.at(2), .site = fields.at(3)};
        searches[search].push_back(
            std::format("{} {} {}", fields.at(5), fields.at(6), fields.at(7)));
    }
    return searches;
}

}  // namespace

// Reference values: PRIMER.exe, True Sites and Almost Sites screens on
// TENZ.GB, TENZC.GB (circular), TLONG.GB and TLONG3.GB
// (tests/legacy/enzymes/find_enzymes.tsv, recorded 2026-10-18)

TEST_CASE("Site search finds the recorded true and almost sites", "[enzymes][legacy]") {
    auto const searches = recorded_searches();
    REQUIRE(searches.size() == 12);
    for (auto const& [search, expected] : searches) {
        auto const molecule = recorded_molecule(search.molecule);
        auto const mismatches = search.kind == "almost" ? Mismatches::at_most_one : Mismatches::none;
        auto const sites = find_sites(molecule.sequence, pattern(search.site), mismatches, molecule.topology);
        auto actual = std::vector<std::string>{};
        std::ranges::transform(sites, std::back_inserter(actual), [&](Site const& site) {
            return describe(site, molecule.sequence, search.site.size(), mismatches);
        });
        INFO(search.molecule << " " << search.kind << " " << search.enzyme);
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
    auto const sites = find_sites(molecule, pattern("GTMKAC"), Mismatches::none, Topology::linear);
    REQUIRE(sites.size() == 2);
    CHECK(sites[0] == Site{.position = 1, .strand = Strand::top, .mismatch_offset = std::nullopt});
    CHECK(sites[1] == Site{.position = 9, .strand = Strand::top, .mismatch_offset = std::nullopt});
}

TEST_CASE("Site search: N matches any base", "[enzymes]") {
    // BglI GCCNNNNNGGC (palindrome)
    auto const molecule = nucleotides("TTGCCATGCAGGCTT");
    auto const sites = find_sites(molecule, pattern("GCCNNNNNGGC"), Mismatches::none, Topology::linear);
    REQUIRE(sites.size() == 1);
    CHECK(sites[0].position == 3);
}

TEST_CASE("Site search: non-palindromic sites on the complement strand", "[enzymes]") {
    // AlwI GGATC; its reverse complement GATCC at 6 is a complement site
    auto const molecule = nucleotides("GGATCGATCC");
    auto const sites = find_sites(molecule, pattern("GGATC"), Mismatches::none, Topology::linear);
    REQUIRE(sites.size() == 2);
    CHECK(sites[0] == Site{.position = 1, .strand = Strand::top, .mismatch_offset = std::nullopt});
    CHECK(sites[1] == Site{.position = 6, .strand = Strand::complement, .mismatch_offset = std::nullopt});
}

TEST_CASE("Site search: one mismatch, and its offset", "[enzymes]") {
    // ApaI GGGCCC: GGGCCA differs at offset 5
    auto const molecule = nucleotides("TTGGGCCATT");
    CHECK(find_sites(molecule, pattern("GGGCCC"), Mismatches::none, Topology::linear).empty());
    auto const sites = find_sites(molecule, pattern("GGGCCC"), Mismatches::at_most_one, Topology::linear);
    REQUIRE(sites.size() == 1);
    CHECK(sites[0] == Site{.position = 3, .strand = Strand::top, .mismatch_offset = 5});
}

TEST_CASE("Site search: a molecule shorter than the site has no site", "[enzymes]") {
    CHECK(find_sites(nucleotides("GAATT"), pattern("GAATTC"), Mismatches::at_most_one, Topology::linear).empty());
    CHECK(find_sites(std::vector<Nucleotide>{}, pattern("GAATTC"), Mismatches::none, Topology::linear).empty());
}

TEST_CASE("Recognition sequences: IUPAC letters only", "[enzymes]") {
    CHECK(parse_site("GAATTC").has_value());
    CHECK(parse_site("gcnngc").has_value());
    CHECK(not parse_site("").has_value());
    CHECK(not parse_site("GAAUUC").has_value());
    CHECK(not parse_site("GAA TTC").has_value());
}

TEST_CASE("Site search: circular molecules, across the origin", "[enzymes]") {
    // EcoRI GAATTC: GAA at the end, TTC at the start of the molecule
    auto const molecule = nucleotides("TTCAAAAAAAGAA");
    CHECK(find_sites(molecule, pattern("GAATTC"), Mismatches::none, Topology::linear).empty());
    auto const sites = find_sites(molecule, pattern("GAATTC"), Mismatches::none, Topology::circular);
    REQUIRE(sites.size() == 1);
    CHECK(sites[0] == Site{.position = 11, .strand = Strand::top, .mismatch_offset = std::nullopt});
}

TEST_CASE("Site search: at the same position, exact sites first", "[enzymes]") {
    // PleI GAGTC; GACTC (its reverse complement) also matches GAGTC with
    // one mismatch: complement exact first, then top with a mismatch
    auto const molecule = nucleotides("TTGACTCTT");
    auto const sites = find_sites(molecule, pattern("GAGTC"), Mismatches::at_most_one, Topology::linear);
    REQUIRE(sites.size() == 2);
    CHECK(sites[0] == Site{.position = 3, .strand = Strand::complement, .mismatch_offset = std::nullopt});
    CHECK(sites[1] == Site{.position = 3, .strand = Strand::top, .mismatch_offset = 2});
}
