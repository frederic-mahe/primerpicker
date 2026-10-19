#include "core/io/pd2_cfg.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/thermo/units.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace {

using Catch::Matchers::WithinAbs;
using primerpicker::io::LegacySystemDefaults;
using primerpicker::io::parse_pd2_cfg;
using primerpicker::io::Pd2CfgError;

auto read_bytes(std::filesystem::path const& path) -> std::vector<std::byte> {
    auto characters = std::vector<char>(std::filesystem::file_size(path));
    auto file = std::ifstream{path, std::ios::binary};
    file.read(characters.data(), std::ssize(characters));
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    auto bytes = std::vector<std::byte>{};
    std::ranges::transform(characters, std::back_inserter(bytes),
                           [](char const character) { return std::bit_cast<std::byte>(character); });
    return bytes;
}

auto settings(char const* const name) -> LegacySystemDefaults {
    auto const result = parse_pd2_cfg(read_bytes(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "config" / name));
    if (not result.defaults.has_value()) {
        throw std::runtime_error{"cannot parse " + std::string{name}};
    }
    return *result.defaults;
}

constexpr auto tolerance_kcal = 1e-6;

}  // namespace

// Reference values: what PRIMER.exe showed in its System Defaults
// forms when it saved these files (tests/legacy/README.md, "Settings
// file pd2.cfg"; tests/legacy/config/, recorded 2026-10-18 and
// 2026-10-19)

TEST_CASE("pd2.cfg: built-in defaults, primer length 25", "[io][legacy]") {
    auto const defaults = settings("PD2_LENGTH25.CFG");
    CHECK(defaults.primer_length == 25);
    // the other values are the built-in defaults (README, "Built-in
    // defaults")
    CHECK(defaults.exclude_tm_at_most_zero);
    CHECK(defaults.sort_homologies_by_tm);
    CHECK(defaults.gc_tolerance_percent == 5);
    CHECK(defaults.tm_tolerance_celsius == 10);
    CHECK(defaults.match_start_offset == 50);
    CHECK(defaults.match_end_offset == 50);
    CHECK(defaults.find_pair_start_offset == 100);
    CHECK(defaults.find_pair_end_offset == 100);
    CHECK(defaults.criteria.gc_min_percent == 50);
    CHECK(defaults.criteria.gc_max_percent == 60);
    CHECK(defaults.criteria.tm_min_celsius == 55);
    CHECK(defaults.criteria.tm_max_celsius == 80);
    CHECK_THAT(defaults.criteria.hairpin_cutoff.kcal_per_mol, WithinAbs(0.0, tolerance_kcal));
    CHECK(defaults.criteria.reject_dimer_3prime_matches == 3);
    CHECK(defaults.criteria.reject_dimer_adjacent_bases == 7);
    CHECK(defaults.criteria.reject_run_length == 3);
    CHECK(defaults.criteria.reject_gc_run_3prime == 3);
    CHECK(defaults.weights.gc_mismatch == 1);
    CHECK(defaults.weights.tm_mismatch == 1);
    CHECK(defaults.weights.dimer_3prime == 1);
    CHECK(defaults.weights.dimer_any == 1);
    CHECK(defaults.weights.false_priming == 1);
}

TEST_CASE("pd2.cfg: changed criteria", "[io][legacy]") {
    // GC 45-65 %, Tm 50-75 °C, hairpin cutoff -1.5 kcal, 3' dimers 4,
    // adjacent bases 8, runs 4, 3' G/C runs 3 (unchanged)
    auto const criteria = settings("PD2_CRITERIA.CFG").criteria;
    CHECK(criteria.gc_min_percent == 45);
    CHECK(criteria.gc_max_percent == 65);
    CHECK(criteria.tm_min_celsius == 50);
    CHECK(criteria.tm_max_celsius == 75);
    CHECK_THAT(criteria.hairpin_cutoff.kcal_per_mol, WithinAbs(-1.5, tolerance_kcal));
    CHECK(criteria.reject_dimer_3prime_matches == 4);
    CHECK(criteria.reject_dimer_adjacent_bases == 8);
    CHECK(criteria.reject_run_length == 4);
    CHECK(criteria.reject_gc_run_3prime == 3);
}

TEST_CASE("pd2.cfg: weights in the order of the form", "[io][legacy]") {
    // typed 1 2 3 4 5, then 3 5 1 4 2, for GC mismatch, Tm mismatch, 3'
    // dimers, any dimers and false priming
    auto const first = settings("PD2_WEIGHTS12345.CFG").weights;
    CHECK(first.gc_mismatch == 1);
    CHECK(first.tm_mismatch == 2);
    CHECK(first.dimer_3prime == 3);
    CHECK(first.dimer_any == 4);
    CHECK(first.false_priming == 5);
    auto const second = settings("PD2_WEIGHTS35142.CFG").weights;
    CHECK(second.gc_mismatch == 3);
    CHECK(second.tm_mismatch == 5);
    CHECK(second.dimer_3prime == 1);
    CHECK(second.dimer_any == 4);
    CHECK(second.false_priming == 2);
}

TEST_CASE("pd2.cfg: a file of another size is rejected", "[io]") {
    auto bytes = read_bytes(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "config" / "PD2_LENGTH25.CFG");
    bytes.pop_back();
    auto const result = parse_pd2_cfg(bytes);
    CHECK(not result.defaults.has_value());
    CHECK(result.error == Pd2CfgError::wrong_size);
}

TEST_CASE("pd2.cfg: a weight above 5 is rejected", "[io]") {
    auto bytes = read_bytes(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "config" / "PD2_LENGTH25.CFG");
    constexpr auto first_weight_offset = 42;
    bytes.at(first_weight_offset) = std::byte{6};
    auto const result = parse_pd2_cfg(bytes);
    CHECK(not result.defaults.has_value());
    CHECK(result.error == Pd2CfgError::invalid_value);
}

TEST_CASE("pd2.cfg: values as written, even unusual ones", "[io]") {
    // config/PD2_RELAXED.CFG (written for the Find Pair probes, with the
    // layout above): GC 0-100 %, Tm 0-99 °C, hairpin cutoff -99 kcal,
    // every reject threshold 99
    auto const criteria = settings("PD2_RELAXED.CFG").criteria;
    CHECK(criteria.gc_min_percent == 0);
    CHECK(criteria.gc_max_percent == 100);
    CHECK(criteria.tm_max_celsius == 99);
    CHECK_THAT(criteria.hairpin_cutoff.kcal_per_mol, WithinAbs(-99.0, tolerance_kcal));
    CHECK(criteria.reject_dimer_3prime_matches == 99);
    CHECK(criteria.reject_gc_run_3prime == 99);
}
