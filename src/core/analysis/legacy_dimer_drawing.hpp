#pragma once

#include "core/analysis/legacy_criteria.hpp"
#include "core/sequence/nucleotide.hpp"

#include <cassert>
#include <cstddef>
#include <span>
#include <string>
#include <string_view>

namespace primerpicker::analysis {

// The duplex of the dimer page of Primer Designer 2.0 (Multi Primer
// Analysis, D; tests/legacy/README.md, "The dimer page"): the alignment
// behind the "3' end matches" value of two primers. The first primer is
// drawn 5' to 3' on the top line, the second 3' to 5' below it, each
// base under its partner, with ':' under the pairs counted from the
// 3' end. The 3' end drawn is the first primer's, unless the second
// primer's count is larger; among alignments of equal count, the one
// facing the largest partner position. Lines start at fixed columns,
// and nothing is drawn left of column 1.
struct DimerDrawing {
    std::string top;
    std::string marks;
    std::string bottom;
};

namespace detail {

constexpr auto dimer_drawing_column = std::size_t{11};  // the first base of the first primer
constexpr auto dimer_five_prime_label = std::string_view{"5'-"};
constexpr auto dimer_three_prime_label = std::string_view{"3'-"};
constexpr auto dimer_three_prime_end = std::string_view{"-3'"};
constexpr auto dimer_five_prime_end = std::string_view{"-5'"};
constexpr auto dimer_first_drawn_column = std::size_t{1};

struct Alignment {
    std::size_t count;             // pairs from the 3' end
    std::size_t partner_position;  // facing the 3' end
};

// The best 3'-end alignment of the first primer on the second (as
// directional_dimer_at_3prime_end), the largest partner position among
// equal counts
constexpr auto best_3prime_alignment(PrimerPair const& pair_of_primers) noexcept -> Alignment {
    auto best = Alignment{.count = 0, .partner_position = 1};
    for (auto partner_position = std::size_t{1}; partner_position < pair_of_primers.second.size(); ++partner_position) {
        auto const count = pairs_from_3prime_end(pair_of_primers, partner_position);
        if (count >= best.count) {
            best = Alignment{.count = count, .partner_position = partner_position};
        }
    }
    return best;
}

// TEXT starting at COLUMN (negative: left of the screen), with nothing
// drawn left of dimer_first_drawn_column
constexpr auto dimer_line(std::string const& text, std::ptrdiff_t const column) -> std::string {
    auto const first = static_cast<std::ptrdiff_t>(dimer_first_drawn_column);
    if (column >= first) {
        return std::string(static_cast<std::size_t>(column), ' ') + text;
    }
    auto const hidden = static_cast<std::size_t>(first - column);
    return std::string(dimer_first_drawn_column, ' ') + (hidden < text.size() ? text.substr(hidden) : std::string{});
}

// Where the duplex goes: the column of the second primer's 3' end, the
// first marked column, and the number of marks
struct Placement {
    std::ptrdiff_t second_3prime_column;
    std::ptrdiff_t marks_column;
    std::size_t count;
};

constexpr auto placement(PrimerPair const& pair_of_primers) noexcept -> Placement {
    auto const& [first, second] = pair_of_primers;
    auto const forward = best_3prime_alignment(pair_of_primers);
    auto const backward = best_3prime_alignment(PrimerPair{.first = second, .second = first});
    auto const column = static_cast<std::ptrdiff_t>(dimer_drawing_column);
    if (backward.count > forward.count) {
        // the second primer's 3' end faces first[partner_position]
        auto const second_3prime_column = column + static_cast<std::ptrdiff_t>(backward.partner_position);
        return {.second_3prime_column = second_3prime_column,
                .marks_column = second_3prime_column,
                .count = backward.count};
    }
    // the first primer's 3' end faces second[partner_position]
    auto const first_3prime_column = column + static_cast<std::ptrdiff_t>(first.size()) - 1;
    auto const second_last = static_cast<std::ptrdiff_t>(second.size()) - 1;
    return {.second_3prime_column =
                first_3prime_column - (second_last - static_cast<std::ptrdiff_t>(forward.partner_position)),
            .marks_column = first_3prime_column - static_cast<std::ptrdiff_t>(forward.count) + 1,
            .count = forward.count};
}

constexpr auto letters(std::span<sequence::Nucleotide const> const primer) -> std::string {
    auto text = std::string{};
    for (auto const nucleotide : primer) {
        text += sequence::to_letter(nucleotide);
    }
    return text;
}

}  // namespace detail

// Both primers need at least two nucleotides. Not noexcept: the strings
// allocate.
[[nodiscard]] constexpr auto legacy_dimer_drawing(PrimerPair const& pair_of_primers) -> DimerDrawing {
    auto const& [first, second] = pair_of_primers;
    assert(first.size() >= 2 and second.size() >= 2);
    auto const column = static_cast<std::ptrdiff_t>(detail::dimer_drawing_column);
    auto const [second_3prime_column, marks_column, count] = detail::placement(pair_of_primers);

    auto drawing = DimerDrawing{};
    drawing.top = detail::dimer_line(
        std::string{detail::dimer_five_prime_label} + detail::letters(first) + std::string{detail::dimer_three_prime_end},
        column - static_cast<std::ptrdiff_t>(detail::dimer_five_prime_label.size()));
    drawing.marks = detail::dimer_line(std::string(count, ':'), marks_column);
    auto reversed = std::string{};
    for (auto index = second.size(); index > 0; --index) {
        reversed += sequence::to_letter(second[index - 1]);
    }
    drawing.bottom = detail::dimer_line(std::string{detail::dimer_three_prime_label} + reversed +
                                            std::string{detail::dimer_five_prime_end},
                                        second_3prime_column -
                                            static_cast<std::ptrdiff_t>(detail::dimer_three_prime_label.size()));
    return drawing;
}

}  // namespace primerpicker::analysis
