#include "core/search/legacy_find_pair.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstddef>
#include <limits>

namespace {

using primerpicker::search::FindPairSettings;
using primerpicker::search::FindPairTargetError;
using primerpicker::search::legacy_check_targets;

struct RecordedCase {
    std::size_t target_a;
    std::size_t target_b;
    std::size_t length;
    FindPairTargetError outcome;
};

constexpr auto test600_size = std::size_t{600};

}  // namespace

// Reference: PRIMER.exe, Find Pair form on pairs/TEST600.GB (600 bases):
// every recorded case of tests/legacy/find_pair_beyond_end/README.txt
// (recorded 2026-10-21; tests/legacy/README.md, "Targets near and
// beyond the molecule's end")
TEST_CASE("Find Pair: the form's checks of the targets", "[search][legacy]") {
    using enum FindPairTargetError;
    constexpr auto recorded = std::array{
        RecordedCase{.target_a = 150, .target_b = 590, .length = 20, .outcome = none},
        RecordedCase{.target_a = 150, .target_b = 600, .length = 20, .outcome = none},
        RecordedCase{.target_a = 150, .target_b = 601, .length = 20, .outcome = target_b_out_of_range},
        RecordedCase{.target_a = 150, .target_b = 650, .length = 20, .outcome = target_b_out_of_range},
        RecordedCase{.target_a = 150, .target_b = 1000, .length = 20, .outcome = target_b_out_of_range},
        RecordedCase{.target_a = 150, .target_b = 600, .length = 20, .outcome = none},
        RecordedCase{.target_a = 579, .target_b = 600, .length = 20, .outcome = none},
        RecordedCase{.target_a = 580, .target_b = 600, .length = 20, .outcome = none},
        RecordedCase{.target_a = 581, .target_b = 600, .length = 20, .outcome = primer_a_too_long},
        RecordedCase{.target_a = 590, .target_b = 600, .length = 20, .outcome = primer_a_too_long},
        RecordedCase{.target_a = 600, .target_b = 600, .length = 20, .outcome = primer_a_too_long},
        RecordedCase{.target_a = 576, .target_b = 600, .length = 25, .outcome = primer_a_too_long},
        RecordedCase{.target_a = 601, .target_b = 600, .length = 20, .outcome = target_a_out_of_range},
        RecordedCase{.target_a = 620, .target_b = 600, .length = 20, .outcome = target_a_out_of_range},
        RecordedCase{.target_a = 650, .target_b = 700, .length = 20, .outcome = target_a_out_of_range},
        RecordedCase{.target_a = 0, .target_b = 600, .length = 20, .outcome = target_a_out_of_range},
        RecordedCase{.target_a = 5, .target_b = 10, .length = 20, .outcome = primer_b_too_long},
        RecordedCase{.target_a = 5, .target_b = 19, .length = 20, .outcome = primer_b_too_long},
        RecordedCase{.target_a = 5, .target_b = 20, .length = 20, .outcome = none},
        RecordedCase{.target_a = 5, .target_b = 21, .length = 20, .outcome = none},
    };
    for (auto const& [target_a, target_b, length, outcome] : recorded) {
        INFO(target_a << ' ' << target_b << ' ' << length);
        auto const settings = FindPairSettings{.target_a = target_a, .target_b = target_b, .length = length};
        CHECK(legacy_check_targets(settings, test600_size) == outcome);
    }
}

// Same reference: the messages of the boxes
TEST_CASE("Find Pair: the messages of the form's checks", "[search][legacy]") {
    using primerpicker::search::legacy_message;
    using enum FindPairTargetError;
    CHECK(legacy_message(target_a_out_of_range) == "Must enter value for this field");
    CHECK(legacy_message(target_b_out_of_range) == "Must enter value for this field");
    CHECK(legacy_message(primer_a_too_long) == "Primer too long to fit within molecule.");
    CHECK(legacy_message(primer_b_too_long) == "Primer too long to fit within molecule.");
}

TEST_CASE("Find Pair: the checks never overflow", "[search][limits]") {
    constexpr auto largest = std::numeric_limits<std::size_t>::max();
    using enum FindPairTargetError;
    CHECK(legacy_check_targets({.target_a = largest, .target_b = 400}, test600_size) == target_a_out_of_range);
    CHECK(legacy_check_targets({.target_a = 16, .target_b = largest}, test600_size) == target_b_out_of_range);
    CHECK(legacy_check_targets({.target_a = 16, .target_b = 400, .length = largest}, test600_size) ==
          primer_a_too_long);
    CHECK(legacy_check_targets({.target_a = 16, .target_b = 400}, largest) == none);
}
