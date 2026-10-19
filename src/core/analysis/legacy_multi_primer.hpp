#pragma once

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_pair.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <span>

namespace primerpicker::analysis {

// The fields Primer A to Primer D of the Multi Primer Analysis of Primer
// Designer 2.0 (tests/legacy/README.md, "Multi Primer Analysis"); a field
// may stay empty
constexpr auto legacy_multi_primer_fields = std::size_t{4};
using MultiPrimerFields = std::array<std::optional<std::span<sequence::Nucleotide const>>, legacy_multi_primer_fields>;

// One column of "Primer Balance and Cautions"
struct LegacyMultiPrimerColumn {
    std::size_t length;
    LegacyActualValues actual;
    LegacyPairCautions cautions;  // the "!" flags, as on the Primer Pair Analysis screen
};

// One cell of the "Primer Dimers" matrices: "3' end" and "Any"
struct LegacyDimerCell {
    std::size_t dimer_3prime_matches;
    std::size_t dimer_adjacent_bases;
    LegacyDimerFlags flags;
};

using LegacyDimerMatrix = std::array<std::array<std::optional<LegacyDimerCell>, legacy_multi_primer_fields>,
                                     legacy_multi_primer_fields>;

struct LegacyMultiPrimer {
    std::array<std::optional<LegacyMultiPrimerColumn>, legacy_multi_primer_fields> columns;
    // dimers[row][column] for row <= column when both fields hold a
    // primer: self-dimers on the diagonal, pairs above it (the cells
    // below it are shown as dots)
    LegacyDimerMatrix dimers;
};

// Multi Primer Analysis: each primer needs at least two nucleotides
[[nodiscard]] constexpr auto legacy_multi_primer(MultiPrimerFields const& fields, LegacyCriteria const& criteria) noexcept
    -> LegacyMultiPrimer {
    auto analysis = LegacyMultiPrimer{};
    for (auto row = std::size_t{0}; row < legacy_multi_primer_fields; ++row) {
        auto const& field = fields.at(row);
        if (not field.has_value()) {
            continue;
        }
        auto const primer = *field;
        auto const actual = legacy_actual_values(primer);
        analysis.columns.at(row) = LegacyMultiPrimerColumn{
            .length = primer.size(), .actual = actual, .cautions = legacy_pair_cautions(actual, criteria)};
        analysis.dimers.at(row).at(row) =
            LegacyDimerCell{.dimer_3prime_matches = actual.dimer_3prime_matches,
                            .dimer_adjacent_bases = actual.dimer_adjacent_bases,
                            .flags = legacy_dimer_flags(actual.dimer_3prime_matches, actual.dimer_adjacent_bases, criteria)};
        for (auto column = row + 1; column < legacy_multi_primer_fields; ++column) {
            auto const& other = fields.at(column);
            if (not other.has_value()) {
                continue;
            }
            auto const values = legacy_pair_values(PrimerPair{.first = primer, .second = *other});
            analysis.dimers.at(row).at(column) = LegacyDimerCell{
                .dimer_3prime_matches = values.dimer_3prime_matches,
                .dimer_adjacent_bases = values.dimer_adjacent_bases,
                .flags = legacy_dimer_flags(values.dimer_3prime_matches, values.dimer_adjacent_bases, criteria)};
        }
    }
    return analysis;
}

}  // namespace primerpicker::analysis
