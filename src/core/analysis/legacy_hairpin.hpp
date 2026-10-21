#pragma once

#include "core/sequence/nucleotide.hpp"
#include "core/thermo/freier_1986.hpp"
#include "core/thermo/units.hpp"

#include <cassert>
#include <cstddef>
#include <string_view>
#include <string>
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
    auto const tenths = hairpin.energy_centi_kcal_per_mol / centi_per_tenth;
    return thermo::KilocaloriesPerMole{tenths / tenths_per_unit};
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

// Whether a Watson-Crick pair further in, reached through Watson-Crick
// and G.T pairs, still closes a loop of at least three bases: the stem
// closed by (inner_5prime, inner_3prime) is then scored there
// (tests/legacy/hairpin_tandem_gt_2/, recorded 2026-10-21)
constexpr auto zips_inward(std::span<sequence::Nucleotide const> const primer, std::size_t const inner_5prime,
                           std::size_t const inner_3prime) noexcept -> bool {
    for (auto depth = std::size_t{1}; inner_5prime + depth < inner_3prime - depth; ++depth) {
        auto const top = primer[inner_5prime + depth];
        auto const bottom = primer[inner_3prime - depth];
        if (inner_3prime - depth - (inner_5prime + depth) - 1 < thermo::freier_1986_min_hairpin_loop) {
            return false;
        }
        if (pairs(top, bottom)) {
            return true;
        }
        if (not is_wobble(top, bottom)) {
            return false;
        }
    }
    return false;
}

}  // namespace detail

// The most stable hairpin of Primer Designer 2.0, if its energy is
// negative (tests/legacy/README.md, "Primer criteria"): stems of
// contiguous pairs, zipped inwards as far as a loop of at
// least three bases allows (over Watson-Crick and G.T pairs, to a
// Watson-Crick pair; recorded 2026-10-21), with G.T
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
            // a stem that could zip further inwards is scored there
            if (detail::zips_inward(primer, inner_5prime, inner_3prime)) {
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

// The drawing of a hairpin on the Primer Analysis screen of Primer
// Designer 2.0 (and on the Multi Primer hairpin page), line by line, as
// on the screen: the 5' part of the primer up to the middle of the loop
// (the larger half of an odd loop), then the pairs (':' under each stem
// pair, ')' one column after the 5' part), then the 3' part, reversed,
// each base under its partner (the middle base of an odd loop is on
// both lines). The 5' part starts at a fixed column; a long 3' part is
// cut at the left edge (nothing is drawn left of column 1);
// tests/legacy/inspect/*.tsv. Not noexcept: the strings allocate.
struct HairpinDrawing {
    std::string top;
    std::string marks;
    std::string bottom;
};

namespace detail {

constexpr auto hairpin_drawing_column = std::size_t{16};  // the first base of the 5' part
constexpr auto five_prime_label = std::string_view{"5'-"};
constexpr auto three_prime_label = std::string_view{"3'-"};
constexpr auto first_drawn_column = std::size_t{1};

// TEXT starting at COLUMN (which may be left of the screen: negative),
// with nothing drawn left of first_drawn_column
constexpr auto clipped_line(std::string const& text, std::ptrdiff_t const column) -> std::string {
    auto const first = static_cast<std::ptrdiff_t>(first_drawn_column);
    if (column >= first) {
        return std::string(static_cast<std::size_t>(column), ' ') + text;
    }
    auto const hidden = static_cast<std::size_t>(first - column);
    return std::string(first_drawn_column, ' ') + (hidden < text.size() ? text.substr(hidden) : std::string{});
}

}  // namespace detail

[[nodiscard]] constexpr auto legacy_hairpin_drawing(std::span<sequence::Nucleotide const> const primer,
                                                    Hairpin const& hairpin) -> HairpinDrawing {
    assert(hairpin.stem_start + (2 * hairpin.stem_length) + hairpin.loop_length <= primer.size());
    auto const column = detail::hairpin_drawing_column;
    auto const top_length = hairpin.stem_start + hairpin.stem_length + ((hairpin.loop_length + 1) / 2);
    // a base of the 3' part shares the column of its partner: i + j is
    // the same for every pair (i, j)
    auto const pair_sum = (2 * hairpin.stem_start) + (2 * hairpin.stem_length) + hairpin.loop_length - 1;
    auto const last = primer.size() - 1;
    // the column of the 3' end (negative when cut at the left edge)
    auto const bottom_column = static_cast<std::ptrdiff_t>(column + pair_sum) - static_cast<std::ptrdiff_t>(last);
    auto const bottom_first = top_length - (hairpin.loop_length % 2);  // the middle base of an odd loop again

    auto drawing = HairpinDrawing{};
    drawing.top = std::string(column - detail::five_prime_label.size(), ' ') + std::string{detail::five_prime_label};
    for (auto index = std::size_t{0}; index < top_length; ++index) {
        drawing.top += sequence::to_letter(primer[index]);
    }
    drawing.marks = std::string(column + hairpin.stem_start, ' ') + std::string(hairpin.stem_length, ':') +
                    std::string(top_length - hairpin.stem_start - hairpin.stem_length, ' ') + ')';
    auto bottom = std::string{detail::three_prime_label};
    for (auto index = last + 1; index > bottom_first; --index) {
        bottom += sequence::to_letter(primer[index - 1]);
    }
    drawing.bottom = detail::clipped_line(
        bottom, bottom_column - static_cast<std::ptrdiff_t>(detail::three_prime_label.size()));
    return drawing;
}

}  // namespace primerpicker::analysis
