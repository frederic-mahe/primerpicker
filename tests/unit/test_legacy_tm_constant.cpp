#include "core/sequence/nucleotide.hpp"
#include "core/thermo/legacy_tm.hpp"

#include <catch2/catch_test_macros.hpp>

#include <charconv>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Recorded {
    std::string primer;
    int tm_celsius;
};

auto number(std::string const& text) -> int {
    auto value = 0;
    auto const* const last = std::next(text.data(), std::ssize(text));
    auto const [end, error] = std::from_chars(text.data(), last, value);
    if (error != std::errc{} or end != last) {
        throw std::runtime_error{"not a number: " + text};
    }
    return value;
}

// tm/tm_constant_probe.tsv: the columns of tm/direct_entry_tm.tsv
auto recordings() -> std::vector<Recorded> {
    constexpr auto column_count = 7U;
    constexpr auto primer_column = 1U;
    constexpr auto tm_column = 5U;
    auto recorded = std::vector<Recorded>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "tm" / "tm_constant_probe.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        if (fields.size() != column_count) {
            throw std::runtime_error{"unexpected line: " + line};
        }
        recorded.push_back({.primer = fields[primer_column], .tm_celsius = number(fields[tm_column])});
    }
    return recorded;
}

}  // namespace

// Sweep (tests/legacy/README.md, "Primer Tm: the constant b"): primers
// whose displayed Tm changes inside the interval left by the earlier
// recordings, typed with Direct Entry
TEST_CASE("Legacy primer Tm: the primers that pin the constant b", "[thermo][legacy]") {
    auto const recorded = recordings();
    REQUIRE(recorded.size() >= 14);
    for (auto const& [primer, tm_celsius] : recorded) {
        INFO(primer);
        auto const sequence = primerpicker::sequence::to_nucleotides(primer);
        REQUIRE(sequence.has_value());
        auto const shown = primerpicker::thermo::legacy_displayed_tm(primerpicker::thermo::legacy_primer_tm(*sequence));
        CHECK(shown == tm_celsius);
    }
}
