#pragma once

#include "core/sequence/genetic_code.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace primerpicker::analysis {

// The translations of the Mutagenize Primer screen of Primer Designer
// 2.0 (tests/legacy/README.md, "Mutagenize"), line by line as inside
// the screen's box: the molecule from 20 bases before the primer to 39
// bases after its 5' end (cut at the molecule's ends), the primer (as
// edited) under its molecule bases, and the reading frames: each amino
// acid under the third base of its codon, in three rows above the
// molecule (O, R, I: codons starting at positions 0, 2 and 1 modulo
// 3, entirely inside the molecule window) and three below the primer
// (M, U, T: positions 1, 2 and 0 modulo 3, entirely inside the primer).
struct MutagenizeView {
    std::array<std::string, 3> original_frames;  // O, R, I
    std::string molecule;
    std::string primer;
    std::array<std::string, 3> mutated_frames;  // M, U, T
};

namespace detail {

constexpr auto mutagenize_primer_column = std::size_t{27};  // the primer's 5' end
constexpr auto mutagenize_bases_before = std::size_t{20};
constexpr auto mutagenize_bases_after = std::size_t{39};  // after the primer's 5' end
constexpr auto codon_length = std::size_t{3};
constexpr auto original_labels = std::array{'O', 'R', 'I'};
constexpr auto original_phases = std::array{std::size_t{0}, std::size_t{2}, std::size_t{1}};
constexpr auto mutated_labels = std::array{'M', 'U', 'T'};
constexpr auto mutated_phases = std::array{std::size_t{1}, std::size_t{2}, std::size_t{0}};

// LETTER at COLUMN of LINE, which grows as needed
constexpr auto put(std::string& line, std::size_t const column, char const letter) -> void {
    if (line.size() <= column) {
        line.resize(column + 1, ' ');
    }
    line[column] = letter;
}

// Where a reading-frame row reads and draws
struct Frame {
    std::size_t first_position;  // 1-based molecule position of the first base read
    std::size_t phase;           // codons start at positions of this value modulo 3
    std::size_t primer_start;    // the position drawn at mutagenize_primer_column
};

// One reading-frame row: the amino acids of the codons of BASES that
// start at a position of the frame's phase, each under its third base
constexpr auto frame_row(char const label, std::span<sequence::Nucleotide const> const bases, Frame const& frame)
    -> std::string {
    auto row = std::string{' '};
    row += label;
    for (auto offset = std::size_t{0}; offset + codon_length <= bases.size(); ++offset) {
        if ((frame.first_position + offset) % codon_length != frame.phase) {
            continue;
        }
        auto const third = frame.first_position + offset + codon_length - 1;  // its position
        put(row, mutagenize_primer_column + third - frame.primer_start,
            sequence::translate(bases[offset], bases[offset + 1], bases[offset + 2]));
    }
    return row;
}

}  // namespace detail

// POSITION: the 1-based molecule position of the primer's 5' end; the
// mutated primer (as long as the original, no 5' extension) must fit
// in the molecule. Not noexcept: the strings allocate.
[[nodiscard]] constexpr auto legacy_mutagenize_view(std::span<sequence::Nucleotide const> const molecule,
                                                    std::size_t const position,
                                                    std::span<sequence::Nucleotide const> const mutated) -> MutagenizeView {
    assert(position >= 1 and position - 1 + mutated.size() <= molecule.size());
    auto const window_first = position > detail::mutagenize_bases_before ? position - detail::mutagenize_bases_before : 1;
    auto const window_last = std::min(molecule.size(), position + detail::mutagenize_bases_after);
    auto const window = molecule.subspan(window_first - 1, window_last - window_first + 1);
    auto const column_of = [position](std::size_t const molecule_position) {
        return detail::mutagenize_primer_column + molecule_position - position;
    };

    auto view = MutagenizeView{};
    for (auto row = std::size_t{0}; row < view.original_frames.size(); ++row) {
        view.original_frames.at(row) = detail::frame_row(
            detail::original_labels.at(row), window,
            {.first_position = window_first, .phase = detail::original_phases.at(row), .primer_start = position});
        view.mutated_frames.at(row) = detail::frame_row(
            detail::mutated_labels.at(row), mutated,
            {.first_position = position, .phase = detail::mutated_phases.at(row), .primer_start = position});
    }
    view.molecule = std::string(column_of(window_first), ' ');
    for (auto const base : window) {
        view.molecule += sequence::to_letter(base);
    }
    view.primer = std::string(detail::mutagenize_primer_column, ' ');
    for (auto const base : mutated) {
        view.primer += sequence::to_letter(base);
    }
    return view;
}

// The strand the Mutagenize Primer screen shows: the primer's own,
// 5' to 3'. For a complement-strand primer, the reverse complement of
// the molecule, on which the primer's 5' end is at molecule size -
// position + 1 (tests/legacy/README.md, "Complement-strand primers").
struct MutagenizeStrand {
    std::vector<sequence::Nucleotide> molecule;
    std::size_t position;  // the primer's 5' end on that strand, 1-based
};

// LOCATION: a primer that lies on the molecule. Not noexcept: the
// vector allocates.
[[nodiscard]] constexpr auto legacy_mutagenize_strand(std::span<sequence::Nucleotide const> const molecule,
                                                      sequence::PrimerLocation const& location) -> MutagenizeStrand {
    assert(location.position >= 1 and location.position <= molecule.size());
    if (location.strand == sequence::Strand::top) {
        return {.molecule = std::vector(molecule.begin(), molecule.end()), .position = location.position};
    }
    auto reversed = std::vector<sequence::Nucleotide>{};
    reversed.reserve(molecule.size());
    std::ranges::transform(molecule.rbegin(), molecule.rend(), std::back_inserter(reversed), sequence::complement);
    return {.molecule = std::move(reversed), .position = molecule.size() - location.position + 1};
}

}  // namespace primerpicker::analysis
