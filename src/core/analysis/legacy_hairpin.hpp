#pragma once

#include "core/sequence/nucleotide.hpp"
#include "core/thermo/freier_1986.hpp"
#include "core/thermo/units.hpp"

#include <cassert>
#include <cstddef>
#include <optional>
#include <span>

namespace primerpicker::analysis {

struct Hairpin {
    std::size_t stem_start;   // 0-based position of the outermost 5' stem base
    std::size_t stem_length;  // base pairs
    std::size_t loop_length;  // unpaired bases between the two arms
    int energy_centi_kcal_per_mol;  // exact: Primer Designer's values are in hundredths
};

// The energy as Primer Designer shows it: truncated toward zero to
// tenths of kcal/mol (a loop of 23 gives -2.16, shown -2.1)
[[nodiscard]] constexpr auto legacy_displayed_energy(Hairpin const& hairpin) noexcept -> thermo::KilocaloriesPerMole {
    constexpr auto centi_per_tenth = 10;
    constexpr auto tenths_per_unit = 10.0;
    // integer division truncates toward zero, as Primer Designer does
    return thermo::KilocaloriesPerMole{(hairpin.energy_centi_kcal_per_mol / centi_per_tenth) / tenths_per_unit};  // NOLINT(bugprone-integer-division)
}

namespace detail {

constexpr auto pairs(sequence::Nucleotide const lhs, sequence::Nucleotide const rhs) noexcept -> bool {
    return sequence::complement(lhs) == rhs;
}

constexpr auto is_wobble(sequence::Nucleotide const lhs, sequence::Nucleotide const rhs) noexcept -> bool {
    using sequence::Nucleotide;
    return (lhs == Nucleotide::guanine and rhs == Nucleotide::thymine) or
           (lhs == Nucleotide::thymine and rhs == Nucleotide::guanine);
}

// Stability (hundredths of kcal/mol) of the step between the pair
// (outer_5prime, outer_3prime) and the next pair inwards; none if the
// step cannot be part of a stem
constexpr auto stack(std::span<sequence::Nucleotide const> const primer, std::size_t const outer_5prime,
                     std::size_t const outer_3prime) noexcept -> std::optional<int> {
    auto const top_5prime = primer[outer_5prime];
    auto const top_3prime = primer[outer_5prime + 1];
    auto const bottom_3prime = primer[outer_3prime];
    auto const bottom_5prime = primer[outer_3prime - 1];
    if (pairs(top_5prime, bottom_3prime) and pairs(top_3prime, bottom_5prime)) {
        return thermo::freier_1986_stack(top_5prime, top_3prime);
    }
    if (not(pairs(top_5prime, bottom_3prime) or is_wobble(top_5prime, bottom_3prime))) {
        return std::nullopt;
    }
    return thermo::freier_1986_wobble_stack(thermo::Step{.top_5prime = thermo::letter(top_5prime),
                                                         .top_3prime = thermo::letter(top_3prime),
                                                         .bottom_3prime = thermo::letter(bottom_3prime),
                                                         .bottom_5prime = thermo::letter(bottom_5prime)});
}

// The hairpin closed by primer[inner_5prime] and primer[inner_3prime],
// with its stem extended outwards as far as bases pair
constexpr auto hairpin_at(std::span<sequence::Nucleotide const> const primer, std::size_t const inner_5prime,
                          std::size_t const inner_3prime) noexcept -> std::optional<Hairpin> {
    auto const loop_length = inner_3prime - inner_5prime - 1;
    auto const loop_penalty = thermo::freier_1986_hairpin_loop(loop_length);
    if (not loop_penalty.has_value()) {
        return std::nullopt;
    }
    // extend over Watson-Crick and G.T pairs (through known steps only),
    // then drop G.T pairs at the outer end: a stem starts and ends with
    // Watson-Crick pairs
    auto length = std::size_t{1};
    while (length <= inner_5prime and inner_3prime + length < primer.size() and
           stack(primer, inner_5prime - length, inner_3prime + length).has_value()) {
        ++length;
    }
    while (length > 1 and not pairs(primer[inner_5prime + 1 - length], primer[inner_3prime + length - 1])) {
        --length;
    }
    auto stem_energy = 0;
    for (auto step = std::size_t{1}; step < length; ++step) {
        stem_energy += stack(primer, inner_5prime - step, inner_3prime + step).value_or(0);
    }
    return Hairpin{.stem_start = inner_5prime + 1 - length,
                   .stem_length = length,
                   .loop_length = loop_length,
                   .energy_centi_kcal_per_mol = *loop_penalty - stem_energy};
}

}  // namespace detail

// The most stable hairpin of Primer Designer 2.0, if its energy is
// negative (tests/legacy/README.md, "Primer criteria"): stems of
// contiguous pairs, zipped inwards as far as a loop of at
// least three bases allows (over Watson-Crick pairs only), with G.T
// pairs inside the stem, scored with Freier 1986 stacks minus the loop
// penalty.
[[nodiscard]] constexpr auto legacy_hairpin(std::span<sequence::Nucleotide const> const primer) noexcept
    -> std::optional<Hairpin> {
    auto best = std::optional<Hairpin>{};
    constexpr auto min_loop = thermo::freier_1986_min_hairpin_loop;
    for (auto inner_5prime = std::size_t{0}; inner_5prime < primer.size(); ++inner_5prime) {
        for (auto inner_3prime = inner_5prime + min_loop + 1; inner_3prime < primer.size(); ++inner_3prime) {
            if (not detail::pairs(primer[inner_5prime], primer[inner_3prime])) {
                continue;
            }
            // a stem that could zip one more pair inwards is scored there
            auto const loop_length = inner_3prime - inner_5prime - 1;
            if (loop_length >= min_loop + 2 and detail::pairs(primer[inner_5prime + 1], primer[inner_3prime - 1])) {
                continue;
            }
            auto const hairpin = detail::hairpin_at(primer, inner_5prime, inner_3prime);
            if (hairpin.has_value() and
                (not best.has_value() or hairpin->energy_centi_kcal_per_mol < best->energy_centi_kcal_per_mol)) {
                best = hairpin;
            }
        }
    }
    if (best.has_value() and best->energy_centi_kcal_per_mol >= 0) {
        return std::nullopt;
    }
    return best;
}

}  // namespace primerpicker::analysis
