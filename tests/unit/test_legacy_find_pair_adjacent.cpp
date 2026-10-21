#include "core/io/molecule_file.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using primerpicker::io::LegacySystemDefaults;
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

auto settings_file(std::filesystem::path const& path) -> LegacySystemDefaults {
    auto const text = read_file(path);
    auto bytes = std::vector<std::byte>{};
    std::ranges::transform(text, std::back_inserter(bytes),
                           [](char const character) { return static_cast<std::byte>(character); });
    auto const result = primerpicker::io::parse_pd2_cfg(bytes);
    if (not result.defaults.has_value()) {
        throw std::runtime_error{"not a pd2.cfg: " + path.string()};
    }
    return *result.defaults;
}

auto trun_molecule() -> std::vector<Nucleotide> {
    auto const text = read_file(legacy_dir / "search" / "TRUN.GB");
    auto const parsed = primerpicker::io::parse_molecule_file({.text = text, .name = "TRUN.GB"});
    auto sequence = parsed.molecule.has_value() ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                                                : std::nullopt;
    if (not sequence.has_value()) {
        throw std::runtime_error{"cannot parse TRUN.GB"};
    }
    return std::move(*sequence);
}

struct Targets {
    std::size_t target_a;
    std::size_t target_b;
};

auto find_pair(std::filesystem::path const& config, Targets const targets) -> primerpicker::search::FindPairResult {
    auto const defaults = settings_file(config);
    auto const settings = primerpicker::search::FindPairSettings{.target_a = targets.target_a,
                                                                 .start_offset_a = defaults.find_pair_start_offset,
                                                                 .target_b = targets.target_b,
                                                                 .end_offset_b = defaults.find_pair_end_offset,
                                                                 .length = defaults.primer_length,
                                                                 .gc_tolerance_percent = defaults.gc_tolerance_percent,
                                                                 .tm_tolerance_celsius = defaults.tm_tolerance_celsius};
    return primerpicker::search::legacy_find_pair(trun_molecule(), settings, defaults.criteria);
}

}  // namespace

// Reference: PRIMER.exe, Find Pair on search/TRUN.GB, targets 1280 and
// 1380, config/PD2_RUNS4.CFG (adjacent-bases criterion 7): 22 and 27
// possible primers, 282 pairs (search/find_pair_3prime_runs.tsv,
// recorded 2026-10-20); 18 of the 300 pairs that meet the other rules
// have 7 A:B adjacent homologous bases (tests/legacy/README.md, "Find
// Pair: A:B adjacent bases", recorded 2026-10-21)
TEST_CASE("Find Pair: fewer A:B adjacent bases than the criterion", "[search][legacy]") {
    auto const result = find_pair(legacy_dir / "config" / "PD2_RUNS4.CFG", {.target_a = 1280, .target_b = 1380});
    CHECK(result.possible_a == 22);
    CHECK(result.possible_b == 27);
    CHECK(result.pairs.size() == 282);
}

// Reference: PRIMER.exe, the same searches with the adjacent-bases
// criterion set to 6, 8 and 5 (find_pair_trun/PD2_RUNS4_ADJ*.CFG,
// inputs; find_pair_trun/README.txt, recorded 2026-10-21; counts
// predicted before the recording)
TEST_CASE("Find Pair: the A:B adjacent-bases rule follows the criterion", "[search][legacy]") {
    auto const directory = legacy_dir / "find_pair_trun";
    auto const adjacent_6 = find_pair(directory / "PD2_RUNS4_ADJ6.CFG", {.target_a = 1280, .target_b = 1380});
    CHECK(adjacent_6.possible_a == 22);
    CHECK(adjacent_6.possible_b == 20);
    CHECK(adjacent_6.pairs.size() == 173);
    auto const adjacent_8 = find_pair(directory / "PD2_RUNS4_ADJ8.CFG", {.target_a = 1280, .target_b = 1380});
    CHECK(adjacent_8.possible_a == 22);
    CHECK(adjacent_8.possible_b == 27);
    CHECK(adjacent_8.pairs.size() == 300);
    auto const adjacent_5 = find_pair(directory / "PD2_RUNS4_ADJ5.CFG", {.target_a = 880, .target_b = 980});
    CHECK(adjacent_5.possible_a == 19);
    CHECK(adjacent_5.possible_b == 29);
    CHECK(adjacent_5.pairs.size() == 160);
}

// Reference: PRIMER.exe, Find Pair on search/TRUN.GB, targets 1680 and
// 1780, config/PD2_RUNS4.CFG: 13 and 24 possible primers, 117 pairs
// (search/find_pair_3prime_runs.tsv, recorded 2026-10-20); primer B at
// 1867 (ATCTGTGGGACCGCCTGCTT) is left out for its hairpin with a tandem
// T.G/G.T step (tests/legacy/README.md, "Find Pair: A:B adjacent bases";
// hairpin_tandem_gt_2/, recorded 2026-10-21)
TEST_CASE("Find Pair: a primer B left out for a hairpin with a tandem G.T step", "[search][legacy]") {
    auto const result = find_pair(legacy_dir / "config" / "PD2_RUNS4.CFG", {.target_a = 1680, .target_b = 1780});
    CHECK(result.possible_a == 13);
    CHECK(result.possible_b == 24);
    CHECK(result.pairs.size() == 117);
}
