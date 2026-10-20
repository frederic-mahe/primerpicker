#pragma once

#include "core/analysis/legacy_builtin_defaults.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/io/pd2.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/thermo/units.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>

namespace primerpicker::io {

// The System Defaults of Primer Designer 2.0, as saved in pd2.cfg
// (tests/legacy/README.md, "Built-in defaults" and "Settings file
// pd2.cfg"); the members' defaults are the built-in ones
struct LegacySystemDefaults {
    std::size_t primer_length = analysis::builtin::primer_length;  // bases, including the 5' extension
    bool exclude_tm_at_most_zero = analysis::builtin::exclude_tm_at_most_zero;  // false priming and uniqueness
    bool sort_homologies_by_tm = analysis::builtin::sort_homologies_by_tm;
    int gc_tolerance_percent = analysis::builtin::gc_tolerance_percent;
    int tm_tolerance_celsius = analysis::builtin::tm_tolerance_celsius;
    std::size_t match_start_offset = analysis::builtin::match_start_offset;  // Match Active Primer: target - offset
    std::size_t match_end_offset = analysis::builtin::match_end_offset;      // to target + offset
    std::size_t find_pair_start_offset = analysis::builtin::find_pair_start_offset;  // Find Pair, primer A: target - offset
    std::size_t find_pair_end_offset = analysis::builtin::find_pair_end_offset;      // primer B: to target + offset
    analysis::LegacyCriteria criteria;
    search::FindPairWeights weights;
};

enum struct Pd2CfgError : std::uint8_t { none, wrong_size, invalid_value };

// C++23 refactoring: replace with std::expected<LegacySystemDefaults, Pd2CfgError>
struct Pd2CfgResult {
    std::optional<LegacySystemDefaults> defaults;
    Pd2CfgError error = Pd2CfgError::none;
};

namespace detail {

constexpr auto pd2_cfg_size = std::size_t{55};
// the first 4 bytes are zero in every recorded file (meaning unknown)
constexpr auto pd2_cfg_unknown_size = std::size_t{4};

// a signed 16-bit value, or none when it is negative
constexpr auto count(ByteReader& reader) noexcept -> std::optional<std::size_t> {
    auto const value = reader.u16();
    if (not value.has_value() or static_cast<std::int16_t>(*value) < 0) {
        return std::nullopt;
    }
    return std::size_t{*value};
}

constexpr auto signed_value(ByteReader& reader) noexcept -> std::optional<int> {
    auto const value = reader.u16();
    return value.has_value() ? std::optional{static_cast<int>(static_cast<std::int16_t>(*value))} : std::nullopt;
}

// a Find Pair weight: 0 (ignored) to 5
constexpr auto max_weight = 5;

constexpr auto weight(ByteReader& reader) noexcept -> std::optional<int> {
    if (reader.remaining() == 0) {
        return std::nullopt;
    }
    auto const value = static_cast<int>(reader.byte_at(0));
    reader.skip(1);
    return value <= max_weight ? std::optional{value} : std::nullopt;
}

}  // namespace detail

// Layout: tests/legacy/README.md, "Settings file pd2.cfg". The weights
// are stored as GC mismatch, Tm mismatch, false priming, 3' dimers, any
// dimers (recorded 2026-10-19); offset 38 holds the run threshold and
// offset 40 the 3' G/C run threshold (config/PD2_RUNS42.CFG, recorded
// 2026-10-19).
[[nodiscard]] constexpr auto parse_pd2_cfg(std::span<std::byte const> const bytes) noexcept -> Pd2CfgResult {
    auto const failure = [](Pd2CfgError const error) { return Pd2CfgResult{.defaults = std::nullopt, .error = error}; };
    if (bytes.size() != detail::pd2_cfg_size) {
        return failure(Pd2CfgError::wrong_size);
    }
    auto reader = detail::ByteReader{bytes};
    reader.skip(detail::pd2_cfg_unknown_size);
    auto const length = detail::count(reader);
    auto const exclude = reader.u16();
    auto const sort = reader.u16();
    auto const gc_tolerance = detail::signed_value(reader);
    auto const tm_tolerance = detail::signed_value(reader);
    // each search range is stored end first, then start (offsets 14 to
    // 20; recorded 2026-10-20, match_variations/PD2_OFFSETS_11_22_33_44)
    auto const match_end = detail::count(reader);
    auto const match_start = detail::count(reader);
    auto const find_pair_end = detail::count(reader);
    auto const find_pair_start = detail::count(reader);
    auto const gc_min = detail::signed_value(reader);
    auto const gc_max = detail::signed_value(reader);
    auto const tm_min = detail::signed_value(reader);
    auto const tm_max = detail::signed_value(reader);
    auto const hairpin = reader.u32();
    auto const dimer_3prime = detail::count(reader);
    auto const dimer_adjacent = detail::count(reader);
    auto const runs = detail::count(reader);
    auto const gc_runs = detail::count(reader);
    auto const gc_weight = detail::weight(reader);
    auto const tm_weight = detail::weight(reader);
    auto const false_priming_weight = detail::weight(reader);
    auto const dimer_3prime_weight = detail::weight(reader);
    auto const dimer_any_weight = detail::weight(reader);
    auto const all_read = length and exclude and sort and gc_tolerance and tm_tolerance and match_start and
                          match_end and find_pair_start and find_pair_end and gc_min and gc_max and tm_min and
                          tm_max and hairpin and dimer_3prime and dimer_adjacent and runs and gc_runs and gc_weight and
                          tm_weight and false_priming_weight and dimer_3prime_weight and dimer_any_weight;
    if (not all_read) {
        return failure(Pd2CfgError::invalid_value);
    }
    return {.defaults =
                LegacySystemDefaults{
                    .primer_length = *length,
                    .exclude_tm_at_most_zero = *exclude != 0,
                    .sort_homologies_by_tm = *sort != 0,
                    .gc_tolerance_percent = *gc_tolerance,
                    .tm_tolerance_celsius = *tm_tolerance,
                    .match_start_offset = *match_start,
                    .match_end_offset = *match_end,
                    .find_pair_start_offset = *find_pair_start,
                    .find_pair_end_offset = *find_pair_end,
                    .criteria = {.gc_min_percent = *gc_min,
                                 .gc_max_percent = *gc_max,
                                 .tm_min_celsius = *tm_min,
                                 .tm_max_celsius = *tm_max,
                                 .hairpin_cutoff = thermo::KilocaloriesPerMole{static_cast<double>(
                                     std::bit_cast<float>(*hairpin))},
                                 .reject_dimer_3prime_matches = *dimer_3prime,
                                 .reject_dimer_adjacent_bases = *dimer_adjacent,
                                 .reject_run_length = *runs,
                                 .reject_gc_run_3prime = *gc_runs},
                    .weights = {.gc_mismatch = *gc_weight,
                                .tm_mismatch = *tm_weight,
                                .dimer_3prime = *dimer_3prime_weight,
                                .dimer_any = *dimer_any_weight,
                                .false_priming = *false_priming_weight},
                },
            .error = Pd2CfgError::none};
}

}  // namespace primerpicker::io
