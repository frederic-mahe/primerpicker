#include "core/io/pd2.hpp"

#include <catch2/catch_test_macros.hpp>

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

using primerpicker::io::parse_pd2;
using primerpicker::io::Pd2Error;
using primerpicker::io::Pd2Primer;
using primerpicker::io::Pd2Result;

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

auto legacy_primer(char const* const name) -> std::vector<std::byte> {
    return read_bytes(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "primers" / name);
}

auto primer_of(Pd2Result const& result) -> Pd2Primer {
    if (not result.primer.has_value()) {
        throw std::runtime_error{"no primer"};
    }
    return *result.primer;
}

}  // namespace

// Reference values: what PRIMER.exe showed in its primer list when it
// saved these files (tests/legacy/primers/, recorded 2026-10-18)

TEST_CASE("a Direct Entry primer with a 5' extension", "[io][legacy]") {
    auto const primer = primer_of(parse_pd2(legacy_primer("M13EXT.PD2")));
    CHECK(primer.extension == "GAATTC");            // PRIMER.exe: GAATTC-GTAAAACGACGGCCAGT
    CHECK(primer.sequence == "GTAAAACGACGGCCAGT");  // PRIMER.exe
    CHECK(primer.molecule_name.empty());            // PRIMER.exe: no molecule
    CHECK_FALSE(primer.complement);
}

TEST_CASE("a primer from a molecule", "[io][legacy]") {
    auto const primer = primer_of(parse_pd2(legacy_primer("COMP.PD2")));
    CHECK(primer.extension.empty());
    CHECK(primer.sequence == "TGTATCCAGGTAATGCAAGG");  // PRIMER.exe
    CHECK(primer.molecule_name == "TEST600");          // PRIMER.exe
    CHECK(primer.molecule_path == "C:\\TEST600.GB");
    CHECK(primer.position == 400);                     // PRIMER.exe: 400(C)
    CHECK(primer.complement);                          // PRIMER.exe: (C)
}

TEST_CASE("saved criteria are those in use", "[io][legacy]") {
    auto const primer = primer_of(parse_pd2(legacy_primer("COMP.PD2")));
    // built-in defaults (tests/legacy/README.md, "Built-in defaults")
    CHECK(primer.criteria.gc_min_percent == 50);
    CHECK(primer.criteria.gc_max_percent == 60);
    CHECK(primer.criteria.tm_min_celsius == 55);
    CHECK(primer.criteria.tm_max_celsius == 80);
    CHECK(primer.criteria.reject_dimer_3prime_matches == 3);
    CHECK(primer.criteria.reject_dimer_adjacent_bases == 7);
    CHECK(primer.criteria.reject_run_length == 3);
    CHECK(primer.criteria.reject_gc_run_3prime == 3);
}

TEST_CASE("a file without the SES header is rejected", "[io]") {
    auto bytes = legacy_primer("COMP.PD2");
    bytes.front() = std::byte{0};
    auto const result = parse_pd2(bytes);
    CHECK_FALSE(result.primer.has_value());
    CHECK(result.error == Pd2Error::bad_header);
}

TEST_CASE("a truncated file is rejected", "[io]") {
    auto bytes = legacy_primer("COMP.PD2");
    constexpr auto truncated_size = 40U;
    bytes.resize(truncated_size);
    auto const result = parse_pd2(bytes);
    CHECK_FALSE(result.primer.has_value());
    CHECK(result.error == Pd2Error::truncated);
}
