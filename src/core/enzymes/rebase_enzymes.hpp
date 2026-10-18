#pragma once

#include "core/enzymes/enzyme.hpp"
#include "core/enzymes/rebase_table.hpp"

#include <algorithm>
#include <functional>
#include <optional>
#include <string_view>

namespace primerpicker::enzymes {

// The enzymes PrimerPicker ships: the commercially available enzymes of
// REBASE (data/rebase/README.md), sorted by name
inline constexpr auto const& rebase_enzymes = rebase_table;

// The enzyme of that name (case sensitive), if any
[[nodiscard]] constexpr auto find_rebase_enzyme(std::string_view const name) noexcept -> std::optional<Enzyme> {
    auto const matches = std::ranges::equal_range(rebase_enzymes, name, std::less{}, &Enzyme::name);
    if (matches.empty()) {
        return std::nullopt;
    }
    return matches.front();
}

}  // namespace primerpicker::enzymes
