#pragma once

#include "core/enzymes/enzyme.hpp"
#include "core/enzymes/rebase_other_table.hpp"
#include "core/enzymes/rebase_table.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iterator>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace primerpicker::enzymes {

// The enzymes PrimerPicker ships: the commercially available enzymes of
// REBASE (data/rebase/README.md), sorted by name
inline constexpr auto const& rebase_enzymes = rebase_table;

// Other REBASE enzymes, accepted as names (old names, enzymes no longer
// sold), sorted by name
inline constexpr auto const& rebase_other_enzymes = rebase_other_table;

enum struct Availability : std::uint8_t { commercial, not_commercial };

struct ResolvedEnzyme {
    Enzyme enzyme;
    Availability availability = Availability::commercial;
};

namespace detail {

constexpr auto find_in(std::span<Enzyme const> const table, std::string_view const name) noexcept
    -> std::optional<Enzyme> {
    auto const matches = std::ranges::equal_range(table, name, std::less{}, &Enzyme::name);
    if (matches.empty()) {
        return std::nullopt;
    }
    return matches.front();
}

struct Abbreviation {
    std::string_view legacy_name;
    std::string_view rebase_name;
};

// Primer Designer 2.0 stores 8 characters per enzyme name
// (tests/legacy/enzymes/enzyme_table.tsv), and shortens these three;
// sorted by legacy name
constexpr auto legacy_abbreviations = std::array{
    Abbreviation{.legacy_name = "Eco47-3", .rebase_name = "Eco47III"},
    Abbreviation{.legacy_name = "HinPI", .rebase_name = "HinP1I"},
    Abbreviation{.legacy_name = "Ksp63", .rebase_name = "Ksp632I"},
};

constexpr auto find_any(std::string_view const name) noexcept -> std::optional<ResolvedEnzyme> {
    if (auto const enzyme = find_in(rebase_enzymes, name); enzyme.has_value()) {
        return ResolvedEnzyme{.enzyme = *enzyme, .availability = Availability::commercial};
    }
    if (auto const enzyme = find_in(rebase_other_enzymes, name); enzyme.has_value()) {
        return ResolvedEnzyme{.enzyme = *enzyme, .availability = Availability::not_commercial};
    }
    return std::nullopt;
}

}  // namespace detail

// The commercially available enzyme of that name (case sensitive), if
// any
[[nodiscard]] constexpr auto find_rebase_enzyme(std::string_view const name) noexcept -> std::optional<Enzyme> {
    return detail::find_in(rebase_enzymes, name);
}

// The REBASE enzyme a name stands for (case sensitive): a commercially
// available enzyme, another REBASE enzyme, one of Primer Designer's
// shortened names, or an old name without its final "I" (Asp718:
// Asp718I, Sau3A: Sau3AI). Not noexcept: the name with "I" allocates.
[[nodiscard]] constexpr auto resolve_enzyme(std::string_view const name) -> std::optional<ResolvedEnzyme> {
    if (name.empty()) {
        return std::nullopt;
    }
    if (auto const found = detail::find_any(name); found.has_value()) {
        return found;
    }
    auto const abbreviation =
        std::ranges::equal_range(detail::legacy_abbreviations, name, std::less{}, &detail::Abbreviation::legacy_name);
    if (not abbreviation.empty()) {
        return detail::find_any(abbreviation.front().rebase_name);
    }
    return detail::find_any(std::string{name} + "I");
}

namespace detail {

template <typename Predicate>
constexpr auto commercial_where(Predicate const predicate) -> std::vector<Enzyme> {
    auto found = std::vector<Enzyme>{};
    std::ranges::copy_if(rebase_enzymes, std::back_inserter(found), predicate);
    return found;
}

}  // namespace detail

// The commercially available enzymes with the same recognition sequence
// and the same cuts (isoschizomers), the enzyme itself included if it
// is sold. Not noexcept: the vector allocates.
[[nodiscard]] constexpr auto commercial_isoschizomers(Enzyme const& enzyme) -> std::vector<Enzyme> {
    return detail::commercial_where([&enzyme](Enzyme const& other) {
        return other.site == enzyme.site and other.cut_count == enzyme.cut_count and other.cuts == enzyme.cuts;
    });
}

// The commercially available enzymes with the same recognition sequence
// and other cuts (neoschizomers). Not noexcept: the vector allocates.
[[nodiscard]] constexpr auto commercial_neoschizomers(Enzyme const& enzyme) -> std::vector<Enzyme> {
    return detail::commercial_where([&enzyme](Enzyme const& other) {
        return other.site == enzyme.site and (other.cut_count != enzyme.cut_count or other.cuts != enzyme.cuts);
    });
}

}  // namespace primerpicker::enzymes
