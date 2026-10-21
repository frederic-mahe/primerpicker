#include "core/glossary/glossary.hpp"
#include "core/glossary/glossary_entry.hpp"
#include "core/glossary/glossary_table.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

using primerpicker::glossary::definition;
using primerpicker::glossary::find_entry;
using primerpicker::glossary::glossary_table;

// every name an entry answers to: its key and its aliases
auto all_names() -> std::multiset<std::string> {
    auto names = std::multiset<std::string>{};
    for (auto const& entry : glossary_table) {
        names.emplace(entry.key);
        auto aliases = entry.aliases;
        while (not aliases.empty()) {
            auto const space = aliases.find(' ');
            names.emplace(aliases.substr(0, space));
            aliases = space == std::string_view::npos ? std::string_view{} : aliases.substr(space + 1);
        }
    }
    return names;
}

// the entry named NAME, which must exist
auto entry_named(std::string_view const name) -> primerpicker::glossary::Entry {
    auto const entry = find_entry(name);
    if (not entry.has_value()) {
        throw std::runtime_error{"no glossary entry " + std::string{name}};
    }
    return *entry;
}

}  // namespace

// Reference: data/glossary/glossary.md (reviewed by the lab, 2026-10-21)
TEST_CASE("the glossary has the 25 reviewed entries", "[glossary]") {
    STATIC_REQUIRE(glossary_table.size() == 25);
    STATIC_REQUIRE(glossary_table.front().key == "primer");
    STATIC_REQUIRE(glossary_table.back().key == "mutagenize");
}

TEST_CASE("every name of the glossary leads to one entry", "[glossary]") {
    auto const names = all_names();
    for (auto const& name : names) {
        INFO(name);
        CHECK(names.count(name) == 1);
        CHECK(std::ranges::none_of(name, [](char const letter) { return letter >= 'A' and letter <= 'Z'; }));
    }
}

TEST_CASE("every entry has a definition", "[glossary]") {
    for (auto const& entry : glossary_table) {
        INFO(entry.key);
        CHECK_FALSE(definition(entry).empty());
        CHECK_FALSE(entry.title.empty());
    }
}

// The decision of 2026-10-20 (proposal §6.5): every criterion and every
// value of the screens has an entry
TEST_CASE("every criterion and screen value has an entry", "[glossary]") {
    constexpr auto required = std::array<std::string_view, 21>{
        "gc",         "tm",           "hairpin",       "dimer",       "runs",     "gc-clamp",
        "criteria",   "extension",    "product-tm",    "annealing-temperature", "balance",
        "uniqueness", "false-priming", "match",        "find-pair",   "restriction-sites",
        "almost-site", "circular",    "multi-primer",  "mutagenize",  "amplicon"};
    for (auto const name : required) {
        INFO(name);
        CHECK(find_entry(name).has_value());
    }
}

TEST_CASE("find_entry: keys and aliases, in any case", "[glossary]") {
    STATIC_REQUIRE(find_entry("GC")->key == "gc");
    STATIC_REQUIRE(find_entry("%GC")->key == "gc");
    STATIC_REQUIRE(find_entry("Ta")->key == "annealing-temperature");
    STATIC_REQUIRE(find_entry("plasmid")->key == "circular");
    STATIC_REQUIRE_FALSE(find_entry("no-such-term").has_value());
    STATIC_REQUIRE_FALSE(find_entry("").has_value());
}

TEST_CASE("definition: the text of the Definition paragraph", "[glossary]") {
    CHECK(definition(entry_named("gc")) ==
          "The percentage of G and C bases in a primer. G·C pairs\n"
          "have three hydrogen bonds and stack more strongly than A·T pairs, so\n"
          "GC-rich primers bind more stably.");
}
