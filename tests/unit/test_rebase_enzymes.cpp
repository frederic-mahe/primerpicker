#include "core/enzymes/enzyme.hpp"
#include "core/enzymes/rebase_enzymes.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

using primerpicker::enzymes::Cut;
using primerpicker::enzymes::Enzyme;
using primerpicker::enzymes::find_rebase_enzyme;
using primerpicker::enzymes::rebase_enzymes;

auto enzyme(std::string_view const name) -> Enzyme {
    auto const found = find_rebase_enzyme(name);
    if (not found.has_value()) {
        throw std::runtime_error{"no enzyme " + std::string{name}};
    }
    return *found;
}

}  // namespace

// Reference values: REBASE version 610, data/rebase/emboss_e.610 and
// emboss_r.610 (the lines are quoted next to each check)

TEST_CASE("REBASE table holds the commercially available enzymes", "[enzymes][rebase]") {
    CHECK(rebase_enzymes.size() == 581);  // emboss_r.610: 581 entries with suppliers
    CHECK(std::ranges::is_sorted(rebase_enzymes, std::less{}, &Enzyme::name));
    CHECK(std::ranges::adjacent_find(rebase_enzymes, std::equal_to{}, &Enzyme::name) == rebase_enzymes.end());
}

TEST_CASE("REBASE enzymes: sites, cuts and suppliers", "[enzymes][rebase]") {
    // emboss_e.610: "EcoRI GAATTC 6 2 0 1 5 0 0"; emboss_r.610 suppliers: BIJNQRSVX
    auto const eco_ri = enzyme("EcoRI");
    CHECK(eco_ri.site == "GAATTC");
    CHECK(eco_ri.cut_count == 2);
    CHECK(eco_ri.cuts[0] == Cut{.top_strand = 1, .bottom_strand = 5});
    CHECK(eco_ri.cuts[1] == Cut{.top_strand = 0, .bottom_strand = 0});
    CHECK(eco_ri.suppliers == "BIJNQRSVX");

    // "MboI GATC 4 2 0 -1 4 0 0": cut before the site
    auto const mbo_i = enzyme("MboI");
    CHECK(mbo_i.cuts[0] == Cut{.top_strand = -1, .bottom_strand = 4});

    // "BglI GCCNNNNNGGC 11 2 0 7 4 0 0": interrupted site
    auto const bgl_i = enzyme("BglI");
    CHECK(bgl_i.site == "GCCNNNNNGGC");

    // "BaeI ACNNNNGTAYC 11 4 0 -11 -16 23 18": two cuts on each strand
    auto const bae_i = enzyme("BaeI");
    CHECK(bae_i.cut_count == 4);
    CHECK(bae_i.cuts[0] == Cut{.top_strand = -11, .bottom_strand = -16});
    CHECK(bae_i.cuts[1] == Cut{.top_strand = 23, .bottom_strand = 18});
}

TEST_CASE("REBASE enzymes: unknown names are not found", "[enzymes][rebase]") {
    CHECK(not find_rebase_enzyme("NoSuchI").has_value());
    CHECK(not find_rebase_enzyme("ecori").has_value());  // names are case sensitive
}
