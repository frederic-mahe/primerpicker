#pragma once

#include "core/glossary/glossary_entry.hpp"
#include "core/glossary/glossary_table.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string_view>

namespace primerpicker::glossary {

namespace detail {

constexpr auto lower(char const letter) noexcept -> char {
    constexpr auto case_offset = 'a' - 'A';
    return letter >= 'A' and letter <= 'Z' ? static_cast<char>(letter + case_offset) : letter;
}

// NAME and CANDIDATE are the same, ASCII case aside
constexpr auto same_name(std::string_view const name, std::string_view const candidate) noexcept -> bool {
    return name.size() == candidate.size() and
           std::ranges::equal(name, candidate, {}, lower, lower);
}

// NAME is one of the space-separated ALIASES
constexpr auto is_alias(std::string_view const name, std::string_view aliases) noexcept -> bool {
    while (not aliases.empty()) {
        auto const space = aliases.find(' ');
        if (same_name(name, aliases.substr(0, space))) {
            return true;
        }
        aliases = space == std::string_view::npos ? std::string_view{} : aliases.substr(space + 1);
    }
    return false;
}

}  // namespace detail

// The entry whose key or alias is NAME (ASCII case aside), if any
[[nodiscard]] constexpr auto find_entry(std::string_view const name) noexcept -> std::optional<Entry> {
    if (name.empty()) {
        return std::nullopt;
    }
    auto const* const found = std::ranges::find_if(glossary_table, [name](Entry const& entry) {
        return detail::same_name(name, entry.key) or detail::is_alias(name, entry.aliases);
    });
    return found == glossary_table.end() ? std::nullopt : std::optional{*found};
}

// The text of the entry's "**Definition.**" paragraph (to its first
// blank line), without that label; empty if there is none
[[nodiscard]] constexpr auto definition(Entry const& entry) noexcept -> std::string_view {
    constexpr auto label = std::string_view{"**Definition.** "};
    auto const start = entry.text.find(label);
    if (start == std::string_view::npos) {
        return {};
    }
    auto const text = entry.text.substr(start + label.size());
    return text.substr(0, text.find("\n\n"));
}

}  // namespace primerpicker::glossary
