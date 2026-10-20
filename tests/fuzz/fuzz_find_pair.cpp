// libFuzzer harness: Find Pair on fuzzed molecules
// (core/search/legacy_find_pair.hpp, and through it the criteria,
// dimer, hairpin and false-priming code). The first bytes give the
// settings and weights, each further byte a base (byte % 4), up to
// max_molecule_length bases.

#include "core/analysis/legacy_verdict.hpp"
#include "core/report/legacy_report.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "fuzz_input.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace {

using primerpicker::search::FindPairSettings;
using primerpicker::search::FindPairWeights;
using primerpicker::search::PairCandidate;
using primerpicker::sequence::Nucleotide;

constexpr auto max_molecule_length = std::size_t{1200};  // the search is slow on long molecules
constexpr auto setting_bytes = std::size_t{14};
constexpr auto max_offset = std::size_t{150};
constexpr auto min_length = std::size_t{8};
constexpr auto length_range = std::size_t{33};  // 8 to 40 bases
constexpr auto tolerance_range = 21;            // 0 to 20
constexpr auto weight_range = 6;                // 0 to 5
constexpr auto byte_bits = 8U;
constexpr auto bases = std::array{Nucleotide::adenine, Nucleotide::cytosine, Nucleotide::guanine, Nucleotide::thymine};

// a position on the molecule, from two bytes (low byte first)
auto position(std::span<std::uint8_t const, 2> const bytes, std::size_t const molecule_size) -> std::size_t {
    auto const value = static_cast<std::size_t>(bytes[0]) | (static_cast<std::size_t>(bytes[1]) << byte_bits);
    return 1 + (value % molecule_size);
}

auto difference(int const lhs, int const rhs) -> int {
    return lhs > rhs ? lhs - rhs : rhs - lhs;
}

}  // namespace

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto const input = primerpicker::fuzz::bytes(data, size);
    if (input.size() <= setting_bytes) {
        return 0;
    }
    auto const settings_input = input.first(setting_bytes);
    auto const molecule_input = input.subspan(setting_bytes).first(
        std::min(input.size() - setting_bytes, max_molecule_length));
    auto molecule = std::vector<Nucleotide>{};
    molecule.reserve(molecule_input.size());
    std::ranges::transform(molecule_input, std::back_inserter(molecule),
                           [](std::uint8_t const byte) { return bases.at(byte % bases.size()); });

    auto const settings = FindPairSettings{
        .target_a = position(settings_input.subspan<0, 2>(), molecule.size()),
        .start_offset_a = settings_input[2] % (max_offset + 1),
        .target_b = position(settings_input.subspan<3, 2>(), molecule.size()),
        .end_offset_b = settings_input[5] % (max_offset + 1),
        .length = min_length + (settings_input[6] % length_range),
        .gc_tolerance_percent = settings_input[7] % tolerance_range,
        .tm_tolerance_celsius = settings_input[8] % tolerance_range};
    auto const weights = FindPairWeights{.gc_mismatch = settings_input[9] % weight_range,
                                         .tm_mismatch = settings_input[10] % weight_range,
                                         .dimer_3prime = settings_input[11] % weight_range,
                                         .dimer_any = settings_input[12] % weight_range,
                                         .false_priming = settings_input[13] % weight_range};

    auto const result = primerpicker::search::legacy_find_pair(molecule, settings, primerpicker::analysis::LegacyCriteria{});
    using primerpicker::fuzz::check;
    check(result.pairs.size() <= result.possible_a * result.possible_b);
    auto const first_a = settings.target_a > settings.start_offset_a ? settings.target_a - settings.start_offset_a : 1;
    for (auto const& pair : result.pairs) {
        check(pair.position_a >= first_a and pair.position_a <= settings.target_a);
        check(pair.position_b >= settings.target_b and pair.position_b <= settings.target_b + settings.end_offset_b);
        check(difference(pair.primer_a.gc_percent, pair.primer_b.gc_percent) <= settings.gc_tolerance_percent);
        check(difference(pair.primer_a.tm_celsius, pair.primer_b.tm_celsius) <= settings.tm_tolerance_celsius);
    }

    // the ranking: same pairs, scores never decrease, ties in generation order
    auto const ranked = primerpicker::search::legacy_rank_pairs(result.pairs, weights);
    check(ranked.size() == result.pairs.size());
    auto const score = [&weights](PairCandidate const& pair) {
        return primerpicker::search::legacy_pair_score(pair, weights);
    };
    for (auto index = std::size_t{1}; index < ranked.size(); ++index) {
        auto const& previous = ranked[index - 1];
        auto const& current = ranked[index];
        check(score(previous) <= score(current));
        if (score(previous) == score(current)) {
            check(previous.position_a < current.position_a or
                  (previous.position_a == current.position_a and previous.position_b < current.position_b));
        }
    }

    // the report: CR LF line ends, one row per printed pair after the
    // list's header
    if (ranked.empty()) {
        return 0;
    }
    auto const best = std::min(ranked.size(), primerpicker::report::legacy_report_best_pairs);
    auto const report = primerpicker::report::legacy_find_pair_report(
        {.molecule_name = "FUZZ", .molecule_size = molecule.size(), .settings = settings, .best_pairs = std::span{ranked}.first(best)},
        primerpicker::analysis::LegacyCriteria{},
        std::chrono::year_month_day{std::chrono::year{2026}, std::chrono::October, std::chrono::day{20}});
    check(primerpicker::fuzz::is_crlf_text(report));
    auto const header_end = report.find("TmA/TmB\r\n\r\n");
    check(header_end != std::string::npos);
    auto const rows = report.substr(header_end + std::string{"TmA/TmB\r\n\r\n"}.size());
    check(std::cmp_equal(std::ranges::count(rows, '\n'), best));
    return 0;
}
