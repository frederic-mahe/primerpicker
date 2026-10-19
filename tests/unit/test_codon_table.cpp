#include "core/sequence/codon_table.hpp"
#include "core/sequence/genetic_code.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>

namespace {

using primerpicker::sequence::legacy_codon_table;
using primerpicker::sequence::Nucleotide;

// "TT(AG),CT(ACGT)" -> TTA, TTG, CTA, CTC, CTG, CTT
auto expand(std::string_view const codons) -> std::set<std::string> {
    auto expanded = std::set<std::string>{};
    auto stream = std::istringstream{std::string{codons}};
    for (auto group = std::string{}; std::getline(stream, group, ',');) {
        auto const open = group.find('(');
        if (open == std::string::npos) {
            expanded.insert(group);
            continue;
        }
        auto const prefix = group.substr(0, open);
        for (auto const third : group.substr(open + 1, group.find(')') - open - 1)) {
            expanded.insert(prefix + third);
        }
    }
    return expanded;
}

constexpr auto bases = std::array{Nucleotide::adenine, Nucleotide::cytosine, Nucleotide::guanine, Nucleotide::thymine};

// the codons of an amino acid in the standard genetic code
auto codons_of(char const code) -> std::set<std::string> {
    auto codons = std::set<std::string>{};
    for (auto const first : bases) {
        for (auto const second : bases) {
            for (auto const third : bases) {
                if (primerpicker::sequence::translate(first, second, third) == code) {
                    codons.insert(std::string{primerpicker::sequence::to_letter(first),
                                              primerpicker::sequence::to_letter(second),
                                              primerpicker::sequence::to_letter(third)});
                }
            }
        }
    }
    return codons;
}

}  // namespace

// Reference values: PRIMER.exe, Mutagenize -> F2
// (tests/legacy/mutagenize/codon_table.tsv, recorded 2026-10-19), and
// the standard genetic code (NCBI translation table 1)

TEST_CASE("the codon table is the recorded one", "[sequence][legacy]") {
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "mutagenize" / "codon_table.tsv"};
    auto recorded = std::vector<std::string>{};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (not line.starts_with('#')) {
            recorded.push_back(line);
        }
    }
    auto table = std::vector<std::string>{};
    std::ranges::transform(legacy_codon_table, std::back_inserter(table), [](auto const& entry) {
        return std::string{entry.name} + '\t' + entry.code + '\t' + std::string{entry.codons};
    });
    CHECK(table == recorded);
}

TEST_CASE("each line of the codon table lists the codons of its amino acid", "[sequence][legacy]") {
    auto all = std::set<std::string>{};
    for (auto const& entry : legacy_codon_table) {
        CAPTURE(entry.name);
        auto const expanded = expand(entry.codons);
        CHECK(expanded == codons_of(entry.code));
        all.insert(expanded.begin(), expanded.end());
    }
    constexpr auto sense_codons = std::size_t{61};
    CHECK(all.size() == sense_codons);
}
