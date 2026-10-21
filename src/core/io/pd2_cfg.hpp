#pragma once

#include "core/analysis/legacy_builtin_defaults.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/io/pd2.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/thermo/units.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace primerpicker::io {

// The bytes of pd2.cfg whose meaning is not known: 0-3 and 47-54
constexpr auto pd2_cfg_header_size = std::size_t{4};
constexpr auto pd2_cfg_trailer_size = std::size_t{8};

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
    // bytes 0-3 and 47-54 of the file: meaning unknown, kept as read when
    // the settings are saved, as PRIMER.exe keeps them (zero when built in)
    std::array<std::byte, pd2_cfg_header_size> unknown_header{};
    std::array<std::byte, pd2_cfg_trailer_size> unknown_trailer{};
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
constexpr auto pd2_cfg_unknown_size = pd2_cfg_header_size;

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
    auto unknown_header = std::array<std::byte, detail::pd2_cfg_unknown_size>{};
    std::ranges::copy(bytes.first(detail::pd2_cfg_unknown_size), unknown_header.begin());
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
    auto unknown_trailer = std::array<std::byte, pd2_cfg_trailer_size>{};
    std::ranges::copy(bytes.last(pd2_cfg_trailer_size), unknown_trailer.begin());
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
                    .unknown_header = unknown_header,
                    .unknown_trailer = unknown_trailer,
                },
            .error = Pd2CfgError::none};
}

namespace detail {

// Little-endian writes of 16-bit values (int16 or count) and bytes
class CfgWriter {
public:
    constexpr void bytes(std::span<std::byte const> const values) {
        bytes_.insert(bytes_.end(), values.begin(), values.end());
    }

    constexpr void u16(std::uint32_t const value) {
        // the constants of pd2.hpp's ByteReader
        bytes_.push_back(static_cast<std::byte>(value & low_byte_mask));
        bytes_.push_back(static_cast<std::byte>((value >> bits_per_byte) & low_byte_mask));
    }

    constexpr void signed_value(int const value) {
        u16(static_cast<std::uint16_t>(static_cast<std::int16_t>(value)));
    }

    constexpr void count(std::size_t const value) { signed_value(static_cast<int>(value)); }

    constexpr void flag(bool const value) { u16(value ? 1U : 0U); }

    constexpr void u32(std::uint32_t const value) {
        u16(value & low_half_mask);
        u16(value >> half_shift);
    }

    constexpr void weight(int const value) { bytes_.push_back(static_cast<std::byte>(value)); }

    [[nodiscard]] constexpr auto result() && -> std::vector<std::byte> { return std::move(bytes_); }

private:
    std::vector<std::byte> bytes_;
};

}  // namespace detail

// The settings file PRIMER.exe writes when System Defaults are saved:
// the layout parse_pd2_cfg reads (55 bytes). Values must fit the file:
// counts and signed values in 16 bits, weights from 0 to 5. Not
// noexcept: the vector allocates.
[[nodiscard]] constexpr auto write_pd2_cfg(LegacySystemDefaults const& defaults) -> std::vector<std::byte> {
    auto writer = detail::CfgWriter{};
    writer.bytes(defaults.unknown_header);
    writer.count(defaults.primer_length);
    writer.flag(defaults.exclude_tm_at_most_zero);
    writer.flag(defaults.sort_homologies_by_tm);
    writer.signed_value(defaults.gc_tolerance_percent);
    writer.signed_value(defaults.tm_tolerance_celsius);
    // each search range is stored end first, then start
    writer.count(defaults.match_end_offset);
    writer.count(defaults.match_start_offset);
    writer.count(defaults.find_pair_end_offset);
    writer.count(defaults.find_pair_start_offset);
    auto const& criteria = defaults.criteria;
    writer.signed_value(criteria.gc_min_percent);
    writer.signed_value(criteria.gc_max_percent);
    writer.signed_value(criteria.tm_min_celsius);
    writer.signed_value(criteria.tm_max_celsius);
    writer.u32(std::bit_cast<std::uint32_t>(static_cast<float>(criteria.hairpin_cutoff.kcal_per_mol)));
    writer.count(criteria.reject_dimer_3prime_matches);
    writer.count(criteria.reject_dimer_adjacent_bases);
    writer.count(criteria.reject_run_length);
    writer.count(criteria.reject_gc_run_3prime);
    // the false-priming weight is stored third
    auto const& weights = defaults.weights;
    writer.weight(weights.gc_mismatch);
    writer.weight(weights.tm_mismatch);
    writer.weight(weights.false_priming);
    writer.weight(weights.dimer_3prime);
    writer.weight(weights.dimer_any);
    writer.bytes(defaults.unknown_trailer);
    return std::move(writer).result();
}

}  // namespace primerpicker::io
