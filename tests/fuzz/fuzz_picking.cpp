// libFuzzer harness: the picking screen and Mutagenize on fuzzed
// molecules (core/search/legacy_picking.hpp, the mutated One Primer
// report of core/report/legacy_report.hpp). The first bytes give the
// strand, the primer's length and position, and the number of keys;
// then come the keys (byte % 8: the six moves and the two directions
// of Find next), then the bases of the molecule (byte % 4), up to
// max_molecule_length. Invariants: the primer always lies on the
// molecule, a refused move keeps it, Find next finds a primer that
// meets the criteria in the right direction, and the screen's lines
// and the report keep their shapes. Added 2026-10-21.

#include "core/analysis/legacy_mutagenize.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/report/legacy_report.hpp"
#include "core/search/legacy_picking.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "fuzz_input.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <span>
#include <string>
#include <vector>

namespace {

using primerpicker::analysis::LegacyCriteria;
using primerpicker::search::PickingKey;
using primerpicker::search::PickingMessage;
using primerpicker::search::SearchDirection;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::PrimerLocation;
using primerpicker::sequence::Strand;
using primerpicker::fuzz::check;

constexpr auto max_molecule_length = std::size_t{2000};  // Find next scans the molecule
constexpr auto header_bytes = std::size_t{5};
constexpr auto max_keys = std::size_t{40};
constexpr auto max_length = std::size_t{60};
constexpr auto byte_bits = 8U;
constexpr auto bases = std::array{Nucleotide::adenine, Nucleotide::cytosine, Nucleotide::guanine, Nucleotide::thymine};
constexpr auto keys = std::array{PickingKey::right, PickingKey::left,   PickingKey::down,
                                 PickingKey::up,    PickingKey::longer, PickingKey::shorter};

auto fits(PrimerLocation const& location, std::span<Nucleotide const> const molecule) -> bool {
    return primerpicker::sequence::primer_at(molecule, location).has_value();
}

auto same(PrimerLocation const& one, PrimerLocation const& other) -> bool {
    return one.position == other.position and one.length == other.length and one.strand == other.strand;
}

auto value_line_shape(std::string const& line) -> bool {
    return line.starts_with("Pos: ") and (line.ends_with("Runs: OK") or line.ends_with("Runs: !"));
}

void check_state(std::span<Nucleotide const> const molecule, PrimerLocation const& location) {
    auto const criteria = LegacyCriteria{};
    auto const state = primerpicker::search::legacy_picking_state(molecule, location, criteria);
    check(state.has_value());
    if (not state.has_value()) {
        return;
    }
    check(value_line_shape(state->value_line));
    // the Mutagenize lines are the value line without its position
    constexpr auto position_columns = std::size_t{12};
    auto const primer = primerpicker::sequence::primer_at(molecule, location);
    check(primer.has_value());
    if (primer.has_value()) {
        check(primerpicker::search::legacy_mutagenize_line(*primer, criteria) ==
              state->value_line.substr(position_columns));
    }
}

void press(std::span<Nucleotide const> const molecule, PrimerLocation& location, std::uint8_t const key_byte) {
    constexpr auto key_count = keys.size() + 2;  // and Find next in both directions
    auto const index = key_byte % key_count;
    if (index < keys.size()) {
        auto const key = keys.at(index);
        auto const [moved, message] = primerpicker::search::legacy_picking_move(location, molecule.size(), key);
        check(fits(moved, molecule));
        // a refused move keeps the primer; only + and - give a message
        if (message != PickingMessage::none) {
            check(same(moved, location));
            check(key == PickingKey::longer or key == PickingKey::shorter);
        }
        location = moved;
        return;
    }
    auto const direction = index == keys.size() ? SearchDirection::forward : SearchDirection::backward;
    auto const found = primerpicker::search::legacy_find_next(molecule, location, LegacyCriteria{}, direction);
    if (not found.has_value()) {
        return;
    }
    check(fits(*found, molecule));
    check(found->length == location.length and found->strand == location.strand);
    check(direction == SearchDirection::forward ? found->position > location.position
                                                : found->position < location.position);
    auto const state = primerpicker::search::legacy_picking_state(molecule, *found, LegacyCriteria{});
    check(state.has_value() and state->meets_criteria);
    location = *found;
}

// Mutagenize on a top-strand primer: the screen, and the mutated report
// of the primer with its bases rotated (a change at most positions)
void check_mutagenize(std::span<Nucleotide const> const molecule, PrimerLocation const& location) {
    if (location.strand != Strand::top or location.length < 2) {
        return;
    }
    auto const original = primerpicker::sequence::primer_at(molecule, location);
    if (not original.has_value()) {
        return;
    }
    auto mutated = *original;
    std::ranges::transform(mutated, mutated.begin(), [](Nucleotide const base) {
        auto const next = (static_cast<std::size_t>(std::ranges::find(bases, base) - bases.begin()) + 1) % bases.size();
        return bases.at(next);
    });
    auto const view = primerpicker::analysis::legacy_mutagenize_view(molecule, location.position, mutated);
    check(not view.molecule.empty() and not view.primer.empty());
    using namespace std::chrono;
    constexpr auto report_year = 2026;
    constexpr auto report_month = 10U;
    constexpr auto report_day = 21U;
    auto const report = primerpicker::report::legacy_mutated_inspect_report(
        {.original = *original,
         .mutated = mutated,
         .molecule_name = "FUZZ",
         .molecule_size = molecule.size(),
         .position = location.position},
        LegacyCriteria{}, year_month_day{year{report_year}, month{report_month}, day{report_day}});
    check(primerpicker::fuzz::is_crlf_text(report));
}

}  // namespace

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto const input = primerpicker::fuzz::bytes(data, size);
    if (input.size() <= header_bytes) {
        return 0;
    }
    auto const strand = (input[0] % 2) == 0 ? Strand::top : Strand::complement;
    auto const length = 1 + (input[1] % max_length);
    auto const position_value = static_cast<std::size_t>(input[2]) | (static_cast<std::size_t>(input[3]) << byte_bits);
    auto const key_count = std::min(static_cast<std::size_t>(input[4]) % (max_keys + 1), input.size() - header_bytes);
    auto const key_bytes = input.subspan(header_bytes, key_count);
    auto const base_bytes = input.subspan(header_bytes + key_count);
    auto molecule = std::vector<Nucleotide>{};
    std::ranges::transform(base_bytes.first(std::min(base_bytes.size(), max_molecule_length)),
                           std::back_inserter(molecule),
                           [](std::uint8_t const byte) { return bases.at(byte % bases.size()); });
    if (molecule.empty()) {
        return 0;
    }
    auto location = PrimerLocation{.position = 1 + (position_value % molecule.size()), .length = length, .strand = strand};
    if (not fits(location, molecule)) {
        return 0;
    }
    check_state(molecule, location);
    for (auto const key_byte : key_bytes) {
        press(molecule, location, key_byte);
        check_state(molecule, location);
    }
    check_mutagenize(molecule, location);
    return 0;
}
