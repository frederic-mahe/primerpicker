// libFuzzer harness: Multi Primer Analysis on fuzzed primers
// (core/analysis/legacy_multi_primer.hpp, and through it the criteria,
// dimer and hairpin code). The first byte says which of the fields A to
// D hold a primer (its low four bits); then, for each filled field, a
// length byte (2 to max_primer_length bases) and the bases (byte % 4).
// Missing bytes leave the remaining fields empty.

#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"
#include "fuzz_input.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <optional>
#include <span>
#include <vector>

namespace {

using primerpicker::analysis::legacy_multi_primer;
using primerpicker::analysis::legacy_multi_primer_fields;
using primerpicker::analysis::LegacyDimerCell;
using primerpicker::analysis::MultiPrimerFields;
using primerpicker::sequence::Nucleotide;

constexpr auto min_primer_length = std::size_t{2};
constexpr auto max_primer_length = std::size_t{60};
constexpr auto bases = std::array{Nucleotide::adenine, Nucleotide::cytosine, Nucleotide::guanine, Nucleotide::thymine};
constexpr auto last_field = legacy_multi_primer_fields - 1;

auto same_cell(std::optional<LegacyDimerCell> const& lhs, std::optional<LegacyDimerCell> const& rhs) -> bool {
    if (not lhs.has_value() or not rhs.has_value()) {
        return lhs.has_value() == rhs.has_value();
    }
    return lhs->dimer_3prime_matches == rhs->dimer_3prime_matches and
           lhs->dimer_adjacent_bases == rhs->dimer_adjacent_bases and
           lhs->flags.three_prime == rhs->flags.three_prime and lhs->flags.adjacent == rhs->flags.adjacent;
}

}  // namespace

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto input = primerpicker::fuzz::bytes(data, size);
    if (input.empty()) {
        return 0;
    }
    auto const filled = input.front();
    input = input.subspan(1);
    auto primers = std::array<std::vector<Nucleotide>, legacy_multi_primer_fields>{};
    auto fields = MultiPrimerFields{};
    for (auto field = std::size_t{0}; field < legacy_multi_primer_fields; ++field) {
        if (((filled >> field) & 1U) == 0 or input.empty()) {
            continue;
        }
        auto const length = min_primer_length + (input.front() % (max_primer_length - min_primer_length + 1));
        input = input.subspan(1);
        if (input.size() < length) {
            break;
        }
        std::ranges::transform(input.first(length), std::back_inserter(primers.at(field)),
                               [](std::uint8_t const byte) { return bases.at(byte % bases.size()); });
        input = input.subspan(length);
        fields.at(field) = primers.at(field);
    }

    auto const criteria = primerpicker::analysis::LegacyCriteria{};
    auto const analysis = legacy_multi_primer(fields, criteria);
    using primerpicker::fuzz::check;
    for (auto row = std::size_t{0}; row < legacy_multi_primer_fields; ++row) {
        auto const& field = fields.at(row);
        check(analysis.columns.at(row).has_value() == field.has_value());
        for (auto column = std::size_t{0}; column < legacy_multi_primer_fields; ++column) {
            // a cell exists on and above the diagonal, for two filled fields
            auto const expected = row <= column and field.has_value() and fields.at(column).has_value();
            check(analysis.dimers.at(row).at(column).has_value() == expected);
        }
        auto const& column = analysis.columns.at(row);
        if (not column.has_value() or not field.has_value()) {
            continue;
        }
        check(column->length == field->size());
        // the diagonal holds the primer's own dimer values
        auto const& self = analysis.dimers.at(row).at(row);
        check(self.has_value() and self->dimer_3prime_matches == column->actual.dimer_3prime_matches and
              self->dimer_adjacent_bases == column->actual.dimer_adjacent_bases);
    }

    // the fields in reverse order give the same cells, mirrored
    auto reversed = MultiPrimerFields{};
    std::ranges::reverse_copy(fields, reversed.begin());
    auto const mirrored = legacy_multi_primer(reversed, criteria);
    for (auto row = std::size_t{0}; row < legacy_multi_primer_fields; ++row) {
        for (auto column = row; column < legacy_multi_primer_fields; ++column) {
            check(same_cell(analysis.dimers.at(row).at(column), mirrored.dimers.at(last_field - column).at(last_field - row)));
        }
    }
    return 0;
}
