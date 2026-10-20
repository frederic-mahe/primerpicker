#include "core/enzymes/site_search.hpp"
#include "core/io/genbank.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "core/sequence/topology.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using primerpicker::enzymes::find_sites;
using primerpicker::enzymes::Mismatches;
using primerpicker::enzymes::site_match_text;
using primerpicker::sequence::Nucleotide;

auto const enzymes_dir = std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "enzymes";

auto read_file(std::filesystem::path const& path) -> std::string {
    auto file = std::ifstream{path, std::ios::binary};
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    auto text = std::ostringstream{};
    text << file.rdbuf();
    return std::move(text).str();
}

struct RecordedMolecule {
    std::vector<Nucleotide> sequence;
    primerpicker::sequence::Topology topology;
};

auto recorded_molecule(std::string const& name) -> RecordedMolecule {
    auto const parsed = primerpicker::io::parse_genbank(read_file(enzymes_dir / name));
    if (not parsed.molecule.has_value()) {
        throw std::runtime_error{"not a GenBank file: " + name};
    }
    auto sequence = primerpicker::sequence::to_nucleotides(parsed.molecule->sequence);
    if (not sequence.has_value()) {
        throw std::runtime_error{"not a DNA molecule: " + name};
    }
    return {.sequence = std::move(*sequence), .topology = parsed.molecule->topology};
}

// A recorded almost site
struct AlmostSite {
    std::string molecule;
    std::string enzyme;
    std::string site;
    std::size_t rank = 0;
    std::size_t position = 0;
    std::string match;
};

auto almost_sites() -> std::vector<AlmostSite> {
    auto table = std::ifstream{enzymes_dir / "find_enzymes.tsv"};
    if (not table) {
        throw std::runtime_error{"cannot read find_enzymes.tsv"};
    }
    auto rows = std::vector<AlmostSite>{};
    for (auto row = std::string{}; std::getline(table, row);) {
        auto fields = std::istringstream{row};
        auto recorded = AlmostSite{};
        auto kind = std::string{};
        auto strand = std::string{};
        fields >> recorded.molecule >> kind >> recorded.enzyme >> recorded.site >> recorded.rank >>
            recorded.position >> strand >> recorded.match;
        if (kind == "almost") {
            rows.push_back(std::move(recorded));
        }
    }
    return rows;
}

auto site_pattern(std::string const& site) -> primerpicker::enzymes::SitePattern {
    auto pattern = primerpicker::enzymes::parse_site(site);
    if (not pattern.has_value()) {
        throw std::runtime_error{"not a recognition sequence: " + site};
    }
    return std::move(*pattern);
}

}  // namespace

// Reference: PRIMER.exe, Molecule-Ops -> Find Enzymes -> Almost Sites:
// the match column of every recorded almost site
// (tests/legacy/enzymes/find_enzymes.tsv, recorded 2026-10-18)
TEST_CASE("site_match_text: the match column of every recorded almost site", "[enzymes][legacy]") {
    auto const recorded = almost_sites();
    for (auto const& almost : recorded) {
        INFO(almost.molecule << ' ' << almost.enzyme << ' ' << almost.rank);
        auto const molecule = recorded_molecule(almost.molecule);
        auto const sites =
            find_sites(molecule.sequence, site_pattern(almost.site), Mismatches::at_most_one, molecule.topology);
        REQUIRE(almost.rank <= sites.size());
        auto const& found = sites.at(almost.rank - 1);
        REQUIRE(found.position == almost.position);
        CHECK(site_match_text(found, molecule.sequence, almost.site.size()) == almost.match);
    }
    CHECK(recorded.size() == 282);
}

TEST_CASE("site_match_text: dots only for a true site", "[enzymes]") {
    auto const molecule = std::vector<Nucleotide>{Nucleotide::adenine, Nucleotide::guanine, Nucleotide::cytosine,
                                                  Nucleotide::thymine};
    auto const site = primerpicker::enzymes::Site{
        .position = 1, .strand = primerpicker::sequence::Strand::top, .mismatch_offset = std::nullopt};
    CHECK(site_match_text(site, molecule, 4) == "....");
}
