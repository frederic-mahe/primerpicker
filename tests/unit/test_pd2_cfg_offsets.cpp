#include "core/io/pd2_cfg.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

auto read_bytes(std::filesystem::path const& path) -> std::vector<std::byte> {
    auto file = std::ifstream{path, std::ios::binary};
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    auto text = std::ostringstream{};
    text << file.rdbuf();
    auto const characters = std::move(text).str();
    auto bytes = std::vector<std::byte>{};
    std::ranges::transform(characters, std::back_inserter(bytes),
                           [](char const character) { return static_cast<std::byte>(character); });
    return bytes;
}

}  // namespace

// match_variations/PD2_OFFSETS_11_22_33_44.CFG holds 11, 22, 33 and 44 at
// offsets 14, 16, 18 and 20; PRIMER.exe's System Defaults screen shows
// "Match Active Primer: Start @ - 22", "End @ + 11", "Find Pair (normal)
// Start @ - 44", "(complement) End @ + 33"
// (match_variations/PD2_OFFSETS_11_22_33_44.screens.txt)
TEST_CASE("pd2.cfg: the search offsets, end before start", "[pd2_cfg][legacy]") {
    auto const result = primerpicker::io::parse_pd2_cfg(read_bytes(
        std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "match_variations" / "PD2_OFFSETS_11_22_33_44.CFG"));
    if (not result.defaults.has_value()) {
        FAIL("cannot parse PD2_OFFSETS_11_22_33_44.CFG");
        return;
    }
    auto const& defaults = *result.defaults;
    CHECK(defaults.match_start_offset == 22);
    CHECK(defaults.match_end_offset == 11);
    CHECK(defaults.find_pair_start_offset == 44);
    CHECK(defaults.find_pair_end_offset == 33);
}
