#include "core/analysis/legacy_verdict.hpp"
#include "core/io/pd2_cfg.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

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

auto criteria(char const* const name) -> primerpicker::analysis::LegacyCriteria {
    auto const result =
        primerpicker::io::parse_pd2_cfg(read_bytes(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "config" / name));
    if (not result.defaults.has_value()) {
        throw std::runtime_error{"cannot parse " + std::string{name}};
    }
    return result.defaults->criteria;
}

}  // namespace

// Reference values: PRIMER.exe saved tests/legacy/config/PD2_RUNS42.CFG
// after "Reject runs of >= 4 bases (any)" and "Reject runs of >= 2 G or
// C at 3'end" were typed in Modify -> System Defaults -> Criteria
// (tests/legacy/README.md, "Settings file pd2.cfg", recorded 2026-10-19)

TEST_CASE("pd2.cfg: the run thresholds, any bases then 3' G/C", "[io][legacy]") {
    auto const runs = criteria("PD2_RUNS42.CFG");
    CHECK(runs.reject_run_length == 4);
    CHECK(runs.reject_gc_run_3prime == 2);
}
