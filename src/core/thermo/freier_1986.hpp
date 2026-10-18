#pragma once

#include "core/sequence/nucleotide.hpp"
#include "core/thermo/units.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <optional>

namespace primerpicker::thermo {

// Freier SM, Kierzek R, Jaeger JA, Sugimoto N, Caruthers MH, Neilson T,
// Turner DH (1986) Improved free-energy parameters for predictions of
// RNA duplex stability. PNAS 83:9373-9377. Free energies at 37 °C, as
// positive stabilities (kcal/mol), read here for DNA (T for U): this
// is what Primer Designer 2.0 uses for hairpins
// (tests/legacy/README.md, "Primer criteria").

// Watson-Crick nearest neighbours, Table 1, indexed by
// [5' nucleotide][3' nucleotide] of a pair on one strand (A, C, G, T)
constexpr auto freier_1986_stack_table = std::array<std::array<double, 4>, 4>{{
    // AA, AC, AG, AT
    {{0.9, 2.1, 1.7, 0.9}},
    // CA, CC, CG, CT
    {{1.8, 2.9, 2.0, 1.7}},
    // GA, GC, GG, GT
    {{2.3, 3.4, 2.9, 2.1}},
    // TA, TC, TG, TT
    {{1.1, 2.3, 1.8, 0.9}},
}};

[[nodiscard]] constexpr auto freier_1986_stack(sequence::Nucleotide const first,
                                               sequence::Nucleotide const second) noexcept
    -> KilocaloriesPerMole {
    auto const row = static_cast<std::size_t>(first);
    auto const column = static_cast<std::size_t>(second);
    assert(row < freier_1986_stack_table.size() and column < freier_1986_stack_table.front().size());
    return KilocaloriesPerMole{freier_1986_stack_table[row][column]};  // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
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
    double stability;
};

constexpr auto freier_1986_wobble_table = std::array{
    WobbleStack{Step{'A', 'G', 'T', 'T'}, 0.5}, WobbleStack{Step{'A', 'T', 'T', 'G'}, 0.7},
    WobbleStack{Step{'C', 'G', 'G', 'T'}, 1.5}, WobbleStack{Step{'C', 'T', 'G', 'G'}, 1.5},
    WobbleStack{Step{'G', 'A', 'T', 'T'}, 0.5}, WobbleStack{Step{'G', 'C', 'T', 'G'}, 1.9},
    WobbleStack{Step{'G', 'G', 'C', 'T'}, 1.3}, WobbleStack{Step{'G', 'G', 'T', 'T'}, 0.5},
    WobbleStack{Step{'G', 'T', 'T', 'G'}, 0.5}, WobbleStack{Step{'T', 'A', 'G', 'T'}, 0.7},
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

// none if the step has no G.T pair (or is not in the table)
[[nodiscard]] constexpr auto freier_1986_wobble_stack(Step const& step) noexcept
    -> std::optional<KilocaloriesPerMole> {
    // NOLINTNEXTLINE(readability-qualified-auto): array iterators are pointers only in some libraries
    auto const found = std::ranges::find_if(freier_1986_wobble_table, [&step](WobbleStack const& entry) {
        return entry.step == step or entry.step == step.reversed();
    });
    if (found == freier_1986_wobble_table.end()) {
        return std::nullopt;
    }
    return KilocaloriesPerMole{found->stability};
}

// Hairpin loop initiation, Table 3, by loop length (penalties, kcal/mol);
// lengths not in the table are interpolated by Primer Designer
// (recorded: 11 -> 4.6, 13 -> 5.3)
constexpr auto freier_1986_min_hairpin_loop = std::size_t{3};
constexpr auto freier_1986_hairpin_loop_table = std::array{
    7.4,  // 3
    5.9,  // 4
    4.4,  // 5
    4.3,  // 6
    4.1,  // 7
    4.1,  // 8
    4.2,  // 9
    4.3,  // 10
    4.6,  // 11, recorded
    4.9,  // 12
    5.3,  // 13, recorded
    5.6,  // 14
};

// none for lengths without a known value
[[nodiscard]] constexpr auto freier_1986_hairpin_loop(std::size_t const length) noexcept
    -> std::optional<KilocaloriesPerMole> {
    if (length < freier_1986_min_hairpin_loop or
        length - freier_1986_min_hairpin_loop >= freier_1986_hairpin_loop_table.size()) {
        return std::nullopt;
    }
    return KilocaloriesPerMole{freier_1986_hairpin_loop_table[length - freier_1986_min_hairpin_loop]};  // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
}

}  // namespace primerpicker::thermo
