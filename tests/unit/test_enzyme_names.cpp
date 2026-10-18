#include "core/enzymes/enzyme.hpp"
#include "core/enzymes/rebase_enzymes.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using primerpicker::enzymes::Availability;
using primerpicker::enzymes::commercial_isoschizomers;
using primerpicker::enzymes::commercial_neoschizomers;
using primerpicker::enzymes::Cut;
using primerpicker::enzymes::Enzyme;
using primerpicker::enzymes::resolve_enzyme;
using primerpicker::enzymes::ResolvedEnzyme;

auto resolved(std::string_view const name) -> ResolvedEnzyme {
    auto const found = resolve_enzyme(name);
    if (not found.has_value()) {
        throw std::runtime_error{"no enzyme " + std::string{name}};
    }
    return *found;
}

auto names(std::vector<Enzyme> const& enzymes) -> std::vector<std::string_view> {
    auto result = std::vector<std::string_view>{};
    std::ranges::transform(enzymes, std::back_inserter(result), &Enzyme::name);
    return result;
}

auto contains(std::vector<std::string_view> const& list, std::string_view const name) -> bool {
    return std::ranges::find(list, name) != list.end();
}

struct LegacyEnzyme {
    std::string name;
    std::string site;
    int cut;
};

// tests/legacy/enzymes/enzyme_table.tsv: the enzymes of PRIMER.exe's
// default menu
auto legacy_menu() -> std::vector<LegacyEnzyme> {
    auto menu = std::vector<LegacyEnzyme>{};
    auto file = std::ifstream{std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "enzymes" / "enzyme_table.tsv"};
    for (auto line = std::string{}; std::getline(file, line);) {
        if (line.starts_with('#')) {
            continue;
        }
        auto stream = std::istringstream{line};
        auto enzyme = LegacyEnzyme{};
        auto overhang = 0;
        auto in_menu = std::string{};
        if (not(stream >> enzyme.name >> enzyme.site >> enzyme.cut >> overhang >> in_menu)) {
            throw std::runtime_error{"cannot read line: " + line};
        }
        if (in_menu == "Y") {
            menu.push_back(enzyme);
        }
    }
    return menu;
}

// PRIMER.exe counts cut positions from 0 (before the site), REBASE's
// EMBOSS files skip 0 (-1 is before the site)
auto rebase_cut(int const legacy_cut) -> int {
    return legacy_cut >= 1 ? legacy_cut : legacy_cut - 1;
}

// the REBASE enzyme of a legacy name has its site and cut
auto check_legacy_enzyme(LegacyEnzyme const& legacy) -> void {
    INFO(legacy.name);
    auto const found = resolve_enzyme(legacy.name);
    REQUIRE(found.has_value());
    auto const& enzyme = found.value_or(ResolvedEnzyme{}).enzyme;
    CHECK(enzyme.site == legacy.site);
    if (enzyme.cut_count != 0) {
        CHECK(enzyme.cuts[0].top_strand == rebase_cut(legacy.cut));
    }
}

}  // namespace

// Reference values: PRIMER.exe's enzyme table
// (tests/legacy/enzymes/enzyme_table.tsv, recorded 2026-10-18)

TEST_CASE("Every enzyme of PRIMER.exe's menu is found in REBASE, with its site and cut", "[enzymes][legacy]") {
    auto const menu = legacy_menu();
    REQUIRE(menu.size() == 89);
    std::ranges::for_each(menu, check_legacy_enzyme);
}

// Reference values: REBASE version 610, data/rebase/emboss_e.610 and
// withrefm.610 (the lines are quoted next to each check)

TEST_CASE("Enzyme names: commercial enzymes", "[enzymes][rebase]") {
    // emboss_e.610: "EcoRI GAATTC 6 2 0 1 5 0 0"
    auto const eco_ri = resolved("EcoRI");
    CHECK(eco_ri.enzyme.name == "EcoRI");
    CHECK(eco_ri.availability == Availability::commercial);
}

TEST_CASE("Enzyme names: old names without the final I", "[enzymes][rebase]") {
    // withrefm.610: "<1>Asp718I", "<3>G^GTACC", no supplier
    auto const asp718 = resolved("Asp718");
    CHECK(asp718.enzyme.name == "Asp718I");
    CHECK(asp718.enzyme.site == "GGTACC");
    CHECK(asp718.enzyme.cuts[0] == Cut{.top_strand = 1, .bottom_strand = 5});
    CHECK(asp718.availability == Availability::not_commercial);
    // emboss_e.610: "Acc65I GGTACC 6 2 0 1 5 0 0", "KpnI GGTACC 6 2 0 5 1 0 0"
    CHECK(contains(names(commercial_isoschizomers(asp718.enzyme)), "Acc65I"));
    CHECK(not contains(names(commercial_isoschizomers(asp718.enzyme)), "KpnI"));
    CHECK(contains(names(commercial_neoschizomers(asp718.enzyme)), "KpnI"));

    // withrefm.610: "<1>Sau3AI", "<3>^GATC"
    auto const sau3a = resolved("Sau3A");
    CHECK(sau3a.enzyme.name == "Sau3AI");
    CHECK(sau3a.enzyme.site == "GATC");
}

TEST_CASE("Enzyme names: REBASE enzymes that are not sold", "[enzymes][rebase]") {
    // withrefm.610: "<1>SciI", "<3>CTC^GAG", no supplier
    auto const sci_i = resolved("SciI");
    CHECK(sci_i.enzyme.site == "CTCGAG");
    CHECK(sci_i.availability == Availability::not_commercial);
    CHECK(sci_i.enzyme.cuts[0] == Cut{.top_strand = 3, .bottom_strand = 3});
}

TEST_CASE("Enzyme names: PRIMER.exe's abbreviations of 8-letter names", "[enzymes][legacy]") {
    // PRIMER.exe names hold 8 characters (enzyme_table.tsv): Eco47-3,
    // HinPI and Ksp63 stand for Eco47III, HinP1I and Ksp632I
    CHECK(resolved("Eco47-3").enzyme.name == "Eco47III");
    CHECK(resolved("HinPI").enzyme.name == "HinP1I");
    CHECK(resolved("Ksp63").enzyme.name == "Ksp632I");
}

TEST_CASE("Enzyme names: unknown names", "[enzymes][rebase]") {
    CHECK(not resolve_enzyme("NoSuchI").has_value());
    CHECK(not resolve_enzyme("ecori").has_value());  // case sensitive
    CHECK(not resolve_enzyme("").has_value());
}
