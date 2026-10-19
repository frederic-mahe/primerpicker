#include "core/analysis/legacy_homology.hpp"
#include "core/enzymes/site_search.hpp"
#include "core/io/molecule_file.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "core/sequence/topology.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using primerpicker::analysis::HomologyKind;
using primerpicker::analysis::legacy_homology_sites;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::Strand;

// PRIMER.exe refuses a sequence text of 64 KB or more ("Sequence too
// large": about 51,000 bases in a GenBank file, 63,000 in a FASTA
// file; tests/legacy/README.md, "Molecule size"); primerpicker has no
// such limit
constexpr auto molecule_length = std::size_t{70'000};
constexpr auto bases_per_line = std::size_t{60};
constexpr auto ecori_position = std::size_t{69'990};  // 1-based
constexpr auto copy_position = std::size_t{68'001};   // 1-based
constexpr auto primer = std::string_view{"CTCAAAGGCTCAGACATCAG"};
constexpr auto primer_3prime_end = std::string_view{"CAGACATCAG"};

// poly-A, with an EcoRI site and a copy of the primer's 3' end planted
auto planted_sequence() -> std::string {
    auto sequence = std::string(molecule_length, 'A');
    sequence.replace(ecori_position - 1, std::string_view{"GAATTC"}.size(), "GAATTC");
    sequence.replace(copy_position - 1, primer_3prime_end.size(), primer_3prime_end);
    return sequence;
}

auto fasta_text(std::string const& sequence) -> std::string {
    auto text = std::string{">LARGE planted sites\n"};
    for (auto start = std::size_t{0}; start < sequence.size(); start += bases_per_line) {
        text += sequence.substr(start, bases_per_line) + '\n';
    }
    return text;
}

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = primerpicker::sequence::to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence"};
    }
    return *converted;
}

auto read_molecule() -> std::vector<Nucleotide> {
    auto const result = primerpicker::io::parse_fasta(fasta_text(planted_sequence()));
    if (not result.molecule.has_value()) {
        throw std::runtime_error{"cannot parse the FASTA text"};
    }
    return nucleotides(result.molecule->sequence);
}

}  // namespace

// Reference values: positions where the sites were planted (constructed
// input)

TEST_CASE("a molecule longer than PRIMER.exe's limit is read whole", "[io][size]") {
    CHECK(read_molecule().size() == molecule_length);
}

TEST_CASE("restriction sites beyond position 65,535 keep their position", "[enzymes][size]") {
    auto const molecule = read_molecule();
    auto const pattern = primerpicker::enzymes::parse_site("GAATTC");
    if (not pattern.has_value()) {
        FAIL("cannot parse GAATTC");
        return;
    }
    auto const sites = primerpicker::enzymes::find_sites(molecule, *pattern, primerpicker::enzymes::Mismatches::none,
                                                         primerpicker::sequence::Topology::linear);
    REQUIRE(not sites.empty());
    CHECK(std::ranges::all_of(sites, [](auto const& site) { return site.position == ecori_position; }));
}

TEST_CASE("false priming sites beyond position 65,535 keep their position", "[analysis][size]") {
    auto const molecule = read_molecule();
    auto const sites = legacy_homology_sites(nucleotides(primer), molecule, HomologyKind::false_priming);
    REQUIRE(sites.size() == 1);
    CHECK(sites.front().position == copy_position);
    CHECK(sites.front().strand == Strand::top);
    CHECK(sites.front().length == primer_3prime_end.size());
}
