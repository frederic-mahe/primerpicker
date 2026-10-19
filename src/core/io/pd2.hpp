#pragma once

#include "core/analysis/legacy_verdict.hpp"
#include "core/thermo/units.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <utility>

namespace primerpicker::io {

// A primer saved by Primer Designer 2.0 (Save Active Primer, .PD2
// file); layout in tests/legacy/README.md, "Primer files .PD2"
struct Pd2Primer {
    std::string extension;      // 5' extension, empty if none
    std::string sequence;       // the primer without its extension
    std::string molecule_name;  // empty for a Direct Entry primer
    std::string molecule_path;
    int position = 0;           // on the molecule (1 for Direct Entry)
    bool complement = false;    // on the complement strand
    analysis::LegacyCriteria criteria;
};

enum struct Pd2Error : std::uint8_t { none, bad_header, truncated, unknown_record };

// C++23 refactoring: replace with std::expected<Pd2Primer, Pd2Error>
struct Pd2Result {
    std::optional<Pd2Primer> primer;
    Pd2Error error = Pd2Error::none;
};

namespace detail {

constexpr auto bits_per_byte = 8U;
constexpr auto low_half_mask = 0xFFFFU;  // the low 16 bits of a 32-bit value
constexpr auto half_shift = 16U;         // the high 16 bits of a 32-bit value
constexpr auto low_byte_mask = 0xFFU;

// Little-endian reads with bounds checks
class ByteReader {
public:
    constexpr explicit ByteReader(std::span<std::byte const> const bytes) noexcept : bytes_{bytes} {}

    [[nodiscard]] constexpr auto remaining() const noexcept -> std::size_t { return bytes_.size() - offset_; }

    constexpr auto skip(std::size_t const count) noexcept -> bool {
        if (remaining() < count) {
            return false;
        }
        offset_ += count;
        return true;
    }

    [[nodiscard]] constexpr auto byte_at(std::size_t const index) const noexcept -> std::uint8_t {
        return std::to_integer<std::uint8_t>(bytes_[offset_ + index]);
    }

    constexpr auto u16() noexcept -> std::optional<std::uint16_t> {
        if (remaining() < 2) {
            return std::nullopt;
        }
        auto const value = static_cast<std::uint16_t>(byte_at(0) | (byte_at(1) << bits_per_byte));
        offset_ += 2;
        return value;
    }

    constexpr auto u32() noexcept -> std::optional<std::uint32_t> {
        constexpr auto size = std::size_t{4};
        if (remaining() < size) {
            return std::nullopt;
        }
        auto value = std::uint32_t{0};
        for (auto index = size; index > 0; --index) {
            value = (value << bits_per_byte) | byte_at(index - 1);
        }
        offset_ += size;
        return value;
    }

    // NUL-terminated text in a field of the given size
    constexpr auto text(std::size_t const size) -> std::optional<std::string> {
        if (remaining() < size) {
            return std::nullopt;
        }
        auto result = std::string{};
        for (auto index = std::size_t{0}; index < size and byte_at(index) != 0; ++index) {
            result.push_back(static_cast<char>(byte_at(index)));
        }
        offset_ += size;
        return result;
    }

private:
    std::span<std::byte const> bytes_;
    std::size_t offset_ = 0;
};

constexpr auto pd2_magic = std::array{std::byte{0x1A}, std::byte{'S'}, std::byte{'E'}, std::byte{'S'}};
constexpr auto pd2_header_size = std::size_t{14};
constexpr auto pd2_end_tag = std::uint16_t{0xFFFF};
constexpr auto text_group = 8U;
constexpr auto number_group = 9U;

enum struct TextField : std::uint8_t { extension = 0, typed_primer = 1, molecule_name = 2, molecule_path = 3, molecule_primer = 4 };
enum struct NumberField : std::uint8_t {
    position = 3,
    complement = 4,
    gc_range = 0x10,
    tm_range = 0x11,
    hairpin_cutoff = 0x12,
    dimer_thresholds = 0x13,
    run_thresholds = 0x14,
};

constexpr auto low_half(std::uint32_t const value) noexcept -> int {
    return static_cast<int>(static_cast<std::int16_t>(value & low_half_mask));
}

constexpr auto high_half(std::uint32_t const value) noexcept -> int {
    return static_cast<int>(static_cast<std::int16_t>(value >> half_shift));
}

constexpr auto set_text(Pd2Primer& primer, TextField const field, std::string value) noexcept -> void {
    switch (field) {
        case TextField::extension:
            primer.extension = std::move(value);
            break;
        case TextField::typed_primer:
        case TextField::molecule_primer:
            primer.sequence = std::move(value);
            break;
        case TextField::molecule_name:
            primer.molecule_name = std::move(value);
            break;
        case TextField::molecule_path:
            primer.molecule_path = std::move(value);
            break;
        default:
            break;  // other text fields are not known yet
    }
}

// Pairs of int16 are stored (second, first) relative to pd2.cfg: GC
// max then min, adjacent bases then 3' matches. The run thresholds
// (both 3 in the recorded files) are read in the same order: 3' G/C
// run, then runs (inferred).
constexpr auto set_number(Pd2Primer& primer, NumberField const field, std::uint32_t const value) noexcept -> void {
    auto& criteria = primer.criteria;
    switch (field) {
        case NumberField::position:
            primer.position = static_cast<int>(value);
            break;
        case NumberField::complement:
            primer.complement = value != 0;
            break;
        case NumberField::gc_range:
            criteria.gc_max_percent = low_half(value);
            criteria.gc_min_percent = high_half(value);
            break;
        case NumberField::tm_range:
            criteria.tm_max_celsius = low_half(value);
            criteria.tm_min_celsius = high_half(value);
            break;
        case NumberField::hairpin_cutoff:
            criteria.hairpin_cutoff = thermo::KilocaloriesPerMole{static_cast<double>(std::bit_cast<float>(value))};
            break;
        case NumberField::dimer_thresholds:
            criteria.reject_dimer_adjacent_bases = static_cast<std::size_t>(low_half(value));
            criteria.reject_dimer_3prime_matches = static_cast<std::size_t>(high_half(value));
            break;
        case NumberField::run_thresholds:
            criteria.reject_gc_run_3prime = static_cast<std::size_t>(low_half(value));
            criteria.reject_run_length = static_cast<std::size_t>(high_half(value));
            break;
        default:
            break;  // lengths (fields 0 to 2) follow from the texts
    }
}

}  // namespace detail

// Not noexcept: the strings allocate
[[nodiscard]] constexpr auto parse_pd2(std::span<std::byte const> const bytes) -> Pd2Result {
    auto const failure = [](Pd2Error const error) { return Pd2Result{.primer = std::nullopt, .error = error}; };
    if (bytes.size() < detail::pd2_header_size or not std::ranges::equal(bytes.first(detail::pd2_magic.size()), detail::pd2_magic)) {
        return failure(Pd2Error::bad_header);
    }
    auto reader = detail::ByteReader{bytes};
    reader.skip(detail::pd2_header_size);
    auto primer = Pd2Primer{};
    while (true) {
        auto const tag = reader.u16();
        if (not tag.has_value()) {
            return failure(Pd2Error::truncated);
        }
        if (*tag == detail::pd2_end_tag) {
            break;
        }
        auto const group = static_cast<unsigned>(*tag >> detail::bits_per_byte);
        auto const field = static_cast<std::uint8_t>(*tag & detail::low_byte_mask);
        if (group == detail::text_group) {
            auto const size = reader.u32();
            auto text = size.has_value() ? reader.text(*size) : std::nullopt;
            if (not text.has_value()) {
                return failure(Pd2Error::truncated);
            }
            detail::set_text(primer, static_cast<detail::TextField>(field), std::move(*text));
            continue;
        }
        if (group == detail::number_group) {
            auto const unused = reader.u32();
            auto const value = reader.u32();
            if (not unused.has_value() or not value.has_value()) {
                return failure(Pd2Error::truncated);
            }
            detail::set_number(primer, static_cast<detail::NumberField>(field), *value);
            continue;
        }
        return failure(Pd2Error::unknown_record);
    }
    return {.primer = std::move(primer), .error = Pd2Error::none};
}

}  // namespace primerpicker::io
