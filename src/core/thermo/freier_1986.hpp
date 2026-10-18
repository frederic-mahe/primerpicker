#pragma once

#include "core/sequence/nucleotide.hpp"
#include "core/thermo/nearest_neighbour_table.hpp"
#include "core/thermo/units.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <optional>

namespace primerpicker::thermo {

// Freier SM, Kierzek R, Jaeger JA, Sugimoto N, Caruthers MH, Neilson T,
// Turner DH (1986) Improved free-energy parameters for predictions of
// RNA duplex stability. PNAS 83:9373-9377. Free energies at 37 °C, as
// positive stabilities (kcal/mol), read here for DNA (T for U): this
// is what Primer Designer 2.0 uses for hairpins
// (tests/legacy/README.md, "Primer criteria").

// Watson-Crick nearest neighbours, Table 1, in hundredths of kcal/mol
// (integers keep sums exact), indexed by [5' nucleotide][3' nucleotide]
// of a pair on one strand (A, C, G, T)
constexpr auto freier_1986_stack_table = NearestNeighbourTable<int>{{
    // AA, AC, AG, AT
    {{90, 210, 170, 90}},
    // CA, CC, CG, CT
    {{180, 290, 200, 170}},
    // GA, GC, GG, GT
    {{230, 340, 290, 210}},
    // TA, TC, TG, TT
    {{110, 230, 180, 90}},
}};

// in hundredths of kcal/mol
[[nodiscard]] constexpr auto freier_1986_stack(sequence::Nucleotide const first,
                                               sequence::Nucleotide const second) noexcept -> int {
    return lookup(freier_1986_stack_table, first, second);
}

// Steps with a G.T pair, as Primer Designer uses them for hairpins. A
// step is written 5'-ab-3' over 3'-cd-5' (a pairs with c, b with d);
// the same step read on the other strand is 5'-dc-3' over 3'-ba-5'.
// These values are fitted on recorded hairpins
// (tests/legacy/README.md): the recordings only determine the sum of
// the two steps around each G.T pair (a G.T pair never ends a stem),
// and these values are one solution, chosen with GG/TT = 0.5.
struct Step {
    char top_5prime;
    char top_3prime;
    char bottom_3prime;  // pairs with top_5prime
    char bottom_5prime;  // pairs with top_3prime

    // the same step, read on the other strand
    [[nodiscard]] constexpr auto reversed() const noexcept -> Step {
        return {.top_5prime = bottom_5prime,
                .top_3prime = bottom_3prime,
                .bottom_3prime = top_3prime,
                .bottom_5prime = top_5prime};
    }

    friend constexpr auto operator==(Step const&, Step const&) noexcept -> bool = default;
};

struct WobbleStack {
    Step step;
    int stability_centi_kcal_per_mol;
};

// The tandem step T.G over G.T is not in the table: stems do not
// extend through it (no recorded stem does)
constexpr auto freier_1986_wobble_table = std::array{
    WobbleStack{Step{'A', 'G', 'T', 'T'}, 50},  WobbleStack{Step{'A', 'T', 'T', 'G'}, 70},
    WobbleStack{Step{'C', 'G', 'G', 'T'}, 150}, WobbleStack{Step{'C', 'T', 'G', 'G'}, 150},
    WobbleStack{Step{'G', 'A', 'T', 'T'}, 50},  WobbleStack{Step{'G', 'C', 'T', 'G'}, 190},
    WobbleStack{Step{'G', 'G', 'C', 'T'}, 130}, WobbleStack{Step{'G', 'G', 'T', 'T'}, 50},
    WobbleStack{Step{'G', 'T', 'T', 'G'}, 50},  WobbleStack{Step{'T', 'A', 'G', 'T'}, 70},
};

constexpr auto letter(sequence::Nucleotide const nucleotide) noexcept -> char {
    switch (nucleotide) {
        case sequence::Nucleotide::adenine:
            return 'A';
        case sequence::Nucleotide::cytosine:
            return 'C';
        case sequence::Nucleotide::guanine:
            return 'G';
        case sequence::Nucleotide::thymine:
            return 'T';
        default:
            return '?';  // unreachable: all enumerators are handled
    }
}

// in hundredths of kcal/mol; none if the step is not in the table
[[nodiscard]] constexpr auto freier_1986_wobble_stack(Step const& step) noexcept -> std::optional<int> {
    // std::array iterators are pointers in libstdc++ and libc++
    auto const* const found = std::ranges::find_if(freier_1986_wobble_table, [&step](WobbleStack const& entry) {
        return entry.step == step or entry.step == step.reversed();
    });
    if (found == freier_1986_wobble_table.end()) {
        return std::nullopt;
    }
    return found->stability_centi_kcal_per_mol;
}

// Hairpin loop initiation, Table 3 (penalties, hundredths of kcal/mol),
// by loop length. Primer Designer interpolates linearly between these
// lengths, and extrapolates beyond 30 with the slope from 20 to 30
// (recorded up to 34; loops of 35 and more never give a stable hairpin
// in a 50-base primer).
struct LoopPenalty {
    std::size_t length;
    int penalty_centi_kcal_per_mol;
};

constexpr auto freier_1986_min_hairpin_loop = std::size_t{3};
constexpr auto freier_1986_hairpin_loop_table = std::array{
    LoopPenalty{3, 740},   LoopPenalty{4, 590},   LoopPenalty{5, 440},   LoopPenalty{6, 430},
    LoopPenalty{7, 410},   LoopPenalty{8, 410},   LoopPenalty{9, 420},   LoopPenalty{10, 430},
    LoopPenalty{12, 490},  LoopPenalty{14, 560},  LoopPenalty{16, 610},  LoopPenalty{18, 670},
    LoopPenalty{20, 710},  LoopPenalty{30, 890},
};

namespace detail {

constexpr auto interpolate(LoopPenalty const& lower, LoopPenalty const& upper, std::size_t const length) noexcept
    -> int {
    auto const span = static_cast<int>(upper.length - lower.length);
    auto const offset = static_cast<int>(length) - static_cast<int>(lower.length);
    // exact for the table: every difference is a multiple of its span
    return lower.penalty_centi_kcal_per_mol +
           ((upper.penalty_centi_kcal_per_mol - lower.penalty_centi_kcal_per_mol) * offset / span);
}

}  // namespace detail

// in hundredths of kcal/mol; none below three bases
[[nodiscard]] constexpr auto freier_1986_hairpin_loop(std::size_t const length) noexcept -> std::optional<int> {
    if (length < freier_1986_min_hairpin_loop) {
        return std::nullopt;
    }
    auto const& table = freier_1986_hairpin_loop_table;
    // the first tabulated length not below this one (the first entry,
    // 3, is the minimum)
    auto const* const upper = std::ranges::find_if(
        std::next(table.cbegin()), table.cend(), [length](LoopPenalty const& entry) { return length <= entry.length; });
    if (upper != table.cend()) {
        return detail::interpolate(*std::prev(upper), *upper, length);
    }
    // beyond the table: extrapolate the last segment
    return detail::interpolate(*std::prev(table.cend(), 2), table.back(), length);
}

}  // namespace primerpicker::thermo
