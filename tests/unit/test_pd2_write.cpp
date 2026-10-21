#include "core/io/pd2.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using primerpicker::io::parse_pd2;
using primerpicker::io::Pd2Primer;
using primerpicker::io::write_pd2;

auto read_bytes(std::filesystem::path const& path) -> std::vector<std::byte> {
    auto file = std::ifstream{path, std::ios::binary};
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    auto text = std::ostringstream{};
    text << file.rdbuf();
    auto const contents = std::move(text).str();
    auto bytes = std::vector<std::byte>{};
    std::ranges::transform(contents, std::back_inserter(bytes),
                           [](char const character) { return static_cast<std::byte>(character); });
    return bytes;
}

auto primer_of(std::vector<std::byte> const& bytes) -> Pd2Primer {
    auto result = parse_pd2(bytes);
    if (not result.primer.has_value()) {
        throw std::runtime_error{"not a .PD2 file"};
    }
    return std::move(*result.primer);
}

}  // namespace

// Reference: the .PD2 files written by PRIMER.exe's Save Active Primer
// (tests/legacy/README.md, "Primer files .PD2", recorded 2026-10-18):
// writing what was read gives the same bytes
TEST_CASE("write_pd2: the recorded files, byte for byte", "[io][legacy]") {
    for (auto const* const name : {"M13EXT.PD2", "COMP.PD2"}) {
        INFO(name);
        auto const recorded = read_bytes(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "primers" / name);
        CHECK(write_pd2(primer_of(recorded)) == recorded);
    }
}

// Same reference: a Direct Entry primer is stored with position 1 and
// without molecule records; texts are NUL-terminated and padded to an
// even length
TEST_CASE("write_pd2: a Direct Entry primer without extension", "[io]") {
    auto primer = Pd2Primer{};
    primer.sequence = "GTAAAACGACGGCCAGTA";  // 18 bases: 19 bytes with the NUL, 20 padded
    primer.position = 1;
    auto const written = write_pd2(primer);
    auto const read = primer_of(written);
    CHECK(read.sequence == primer.sequence);
    CHECK(read.extension.empty());
    CHECK(read.molecule_name.empty());
    CHECK(read.position == 1);
    CHECK_FALSE(read.complement);
    CHECK(written.size() % 2 == 0);
}
