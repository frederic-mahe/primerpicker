#pragma once

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_hairpin.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>

// The picking screen of Primer Designer 2.0 (MODIFY PRIMER POSITION,
// opened by Primer/Create -> From Molecule and Modify -> Position): the
// primer is a window on the molecule, moved and resized with keys, its
// values shown on one line and updated on every move
// (tests/legacy/README.md, "The picking screen").

namespace primerpicker::search {

// The keys that move or resize the window: the arrows (one base), Up
// and Down (one line of the screen, 60 bases), + and - (one base more
// or less). On both strands, Right and Down raise the position.
enum struct PickingKey : std::uint8_t { right, left, down, up, longer, shorter };

enum struct PickingMessage : std::uint8_t {
    none,
    primer_too_long,          // "Primer Too Long - Extends outside molecule"
    primer_must_have_length,  // "Primer must have length!" (- on one base)
};

struct PickingMove {
    sequence::PrimerLocation location;
    PickingMessage message;
};

namespace detail {

constexpr auto picking_line_bases = std::size_t{60};

// whether the primer lies on a molecule of molecule_size bases
constexpr auto fits(sequence::PrimerLocation const& location, std::size_t const molecule_size) noexcept -> bool {
    auto const& [position, length, strand] = location;
    if (position == 0 or length == 0 or position > molecule_size) {
        return false;
    }
    if (strand == sequence::Strand::top) {
        return length <= molecule_size - (position - 1);
    }
    return length <= position;
}

constexpr auto append_number(std::string& text, std::size_t value) -> void {
    constexpr auto base = std::size_t{10};
    if (value == 0) {
        text += '0';
        return;
    }
    auto digits = std::string{};
    for (; value > 0; value /= base) {
        digits.insert(digits.begin(), static_cast<char>('0' + (value % base)));
    }
    text += digits;
}

// a number right-aligned in 3 columns, as the %GC and the Tm
constexpr auto append_value(std::string& text, std::size_t const value) -> void {
    constexpr auto width = std::size_t{3};
    auto number = std::string{};
    append_number(number, value);
    if (number.size() < width) {
        text.append(width - number.size(), ' ');
    }
    text += number;
}

constexpr auto pad_to(std::string& text, std::size_t const width) -> void {
    if (text.size() < width) {
        text.append(width - text.size(), ' ');
    }
}

// The values of a one-base primer, which the picking screen allows:
// those of legacy_actual_values, with an unknown Tm
[[nodiscard]] constexpr auto single_base_values(std::span<sequence::Nucleotide const> const primer)
    -> analysis::LegacyActualValues {
    assert(primer.size() == 1);
    auto const hairpin = analysis::legacy_hairpin(primer);
    return {.gc_percent = analysis::legacy_gc_percent(primer),
            .tm_celsius = 0,
            .hairpin = hairpin.has_value() ? std::optional{analysis::legacy_displayed_energy(*hairpin)} : std::nullopt,
            .dimer_3prime_matches = analysis::self_dimer_at_3prime_end(primer),
            .dimer_adjacent_bases = analysis::longest_self_complementary_run(primer),
            .longest_run = analysis::longest_run_before_3prime_end(primer),
            .gc_run_3prime = analysis::gc_run_at_3prime_end(primer)};
}

// The values of the picking screen for a primer of one base or more.
// Unlike the Primer Analysis screen, the run at the 3' end counts as a
// run, as in Find Pair (tests/legacy/README.md, "The picking screen").
[[nodiscard]] constexpr auto picking_values(std::span<sequence::Nucleotide const> const primer)
    -> analysis::LegacyActualValues {
    auto values = primer.size() == 1 ? single_base_values(primer) : analysis::legacy_actual_values(primer);
    values.longest_run = std::max(values.longest_run, analysis::run_at_3prime_end(primer));
    return values;
}

}  // namespace detail

// The window after a key: a move that would take the primer off the
// molecule is ignored (no message); a longer primer that would not fit
// is refused with "Primer Too Long", a primer of no base with "Primer
// must have length!". Requires a primer that fits.
[[nodiscard]] constexpr auto legacy_picking_move(sequence::PrimerLocation const& location,
                                                 std::size_t const molecule_size, PickingKey const key) noexcept
    -> PickingMove {
    assert(detail::fits(location, molecule_size));
    auto moved = location;
    switch (key) {
        case PickingKey::right:
            ++moved.position;
            break;
        case PickingKey::left:
            --moved.position;  // 0 does not fit
            break;
        case PickingKey::down:
            // beyond the molecule: a position that does not fit (no overflow)
            moved.position = molecule_size - moved.position >= detail::picking_line_bases
                                 ? moved.position + detail::picking_line_bases
                                 : 0;
            break;
        case PickingKey::up:
            // below the first line: a position that does not fit
            moved.position = moved.position > detail::picking_line_bases
                                 ? moved.position - detail::picking_line_bases
                                 : 0;
            break;
        case PickingKey::longer:
            ++moved.length;
            break;
        case PickingKey::shorter:
            --moved.length;
            break;
        default:
            break;  // unreachable: all enumerators are handled
    }
    if (detail::fits(moved, molecule_size)) {
        return {.location = moved, .message = PickingMessage::none};
    }
    auto const message = [key] {
        switch (key) {
            case PickingKey::longer:
                return PickingMessage::primer_too_long;
            case PickingKey::shorter:
                return PickingMessage::primer_must_have_length;
            case PickingKey::right:
            case PickingKey::left:
            case PickingKey::down:
            case PickingKey::up:
            default:
                return PickingMessage::none;
        }
    }();
    return {.location = location, .message = message};
}

// Find next: Ctrl+Right searches toward higher positions, Ctrl+Left
// toward lower ones, on the same strand
enum struct SearchDirection : std::uint8_t { forward, backward };

// The nearest primer of the same length and strand, beyond the given
// one in the direction, that meets every criterion ("MEETS CRITERIA
// SELECTED: Yes"); none: PRIMER.exe's "No Primers found" (the window
// stays). Not noexcept: primers are copied into vectors.
[[nodiscard]] constexpr auto legacy_find_next(std::span<sequence::Nucleotide const> const molecule,
                                              sequence::PrimerLocation const& from,
                                              analysis::LegacyCriteria const& criteria,
                                              SearchDirection const direction)
    -> std::optional<sequence::PrimerLocation> {
    auto candidate = from;
    while (true) {
        if (direction == SearchDirection::forward) {
            if (candidate.position >= molecule.size()) {
                return std::nullopt;
            }
            ++candidate.position;
        } else {
            if (candidate.position <= 1) {
                return std::nullopt;
            }
            --candidate.position;
        }
        auto const primer = sequence::primer_at(molecule, candidate);
        if (not primer.has_value()) {
            continue;
        }
        if (analysis::evaluate(detail::picking_values(*primer), criteria).all()) {
            return candidate;
        }
    }
}

// What the picking screen shows for a primer
struct PickingState {
    // "Pos: 51     %GC: 50    Tm°C: 66    Hairpin: OK  Dimers: OK  Runs: !"
    std::string value_line;
    // "MEETS CRITERIA SELECTED: Yes"; the window is then green (light
    // grey otherwise)
    bool meets_criteria = false;
};


// The value line (tests/legacy/README.md, "The picking screen"): "Pos: "
// and the position in 7 columns; "%GC:" and "Tm°C:" with the value
// right-aligned in 3 columns ("?" in the middle one when the Tm is not
// positive), then '!' when out of range or a space, then 3 spaces; then
// OK or ! (in 2 columns, then 2 spaces) for the hairpin, the dimers (3'
// end and adjacent bases) and the runs (any base, the run at the 3' end
// included, and G or C at the 3' end). Not noexcept: builds a string.
[[nodiscard]] constexpr auto legacy_value_line(std::size_t const position,
                                               analysis::LegacyActualValues const& actual,
                                               analysis::LegacyVerdict const& verdict) -> std::string {
    using namespace std::string_view_literals;
    constexpr auto position_end = std::size_t{12};  // "Pos: " and 7 columns
    constexpr auto spaces = "   "sv;
    constexpr auto word_width = std::size_t{4};
    auto const flag = [](bool const meets) { return meets ? ' ' : '!'; };
    auto const word = [](std::string_view const label, bool const meets) {
        auto field = std::string{label};
        field += meets ? "OK"sv : "!"sv;
        detail::pad_to(field, label.size() + word_width);
        return field;
    };
    auto line = std::string{"Pos: "};
    detail::append_number(line, position);
    detail::pad_to(line, position_end);
    line += "%GC:";
    detail::append_value(line, static_cast<std::size_t>(actual.gc_percent));
    line += flag(verdict.gc);
    line += spaces;
    line += "Tm°C:";
    if (actual.tm_celsius > 0) {
        detail::append_value(line, static_cast<std::size_t>(actual.tm_celsius));
    } else {
        line += " ? ";
    }
    line += flag(verdict.tm);
    line += spaces;
    line += word("Hairpin: ", verdict.no_hairpin);
    line += word("Dimers: ", verdict.no_3prime_dimer and verdict.no_dimer);
    line += "Runs: ";
    line += verdict.no_runs and verdict.no_3prime_gc_run ? "OK"sv : "!"sv;
    return line;
}

// The picking screen's values for the primer at a location; none if it
// does not fit on the molecule. Not noexcept: the primer is copied into
// a vector.
[[nodiscard]] constexpr auto legacy_picking_state(std::span<sequence::Nucleotide const> const molecule,
                                                  sequence::PrimerLocation const& location,
                                                  analysis::LegacyCriteria const& criteria)
    -> std::optional<PickingState> {
    auto const primer = sequence::primer_at(molecule, location);
    if (not primer.has_value()) {
        return std::nullopt;
    }
    auto const actual = detail::picking_values(*primer);
    auto const verdict = analysis::evaluate(actual, criteria);
    return PickingState{.value_line = legacy_value_line(location.position, actual, verdict),
                        .meets_criteria = verdict.all()};
}

}  // namespace primerpicker::search
