#pragma once

#include "core/io/genbank.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace primerpicker::io {

// Molecule files that Primer Designer 2.0 loads (tests/legacy/README.md,
// "Molecule file formats"): GenBank, FASTA (name: the first word of the
// header) and plain sequence (name: the file name). It rejects EMBL.

enum struct MoleculeError : std::uint8_t { none, missing_header, missing_origin, empty_sequence };

// C++23 refactoring: replace with std::expected<Molecule, MoleculeError>
struct MoleculeResult {
    std::optional<Molecule> molecule;
    MoleculeError error = MoleculeError::none;
};

// A file's contents and its name (a plain sequence is named after its
// file)
struct MoleculeFile {
    std::string_view text;
    std::string_view name;
};

namespace detail {

constexpr auto failure(MoleculeError const error) -> MoleculeResult {
    return {.molecule = std::nullopt, .error = error};
}

constexpr auto success(std::string name, std::string sequence) -> MoleculeResult {
    if (sequence.empty()) {
        return failure(MoleculeError::empty_sequence);
    }
    return {.molecule = Molecule{.name = std::move(name), .sequence = std::move(sequence)}, .error = MoleculeError::none};
}

constexpr auto is_blank_line(std::string_view const line) noexcept -> bool {
    return line.find_first_not_of(" \t") == std::string_view::npos;
}

// the first line that is not blank, if any
constexpr auto first_content_line(std::string_view const text) noexcept -> std::optional<std::string_view> {
    auto lines = Lines{text};
    for (auto line = lines.next(); line.has_value(); line = lines.next()) {
        if (not is_blank_line(*line)) {
            return line;
        }
    }
    return std::nullopt;
}

}  // namespace detail

// The first record of a FASTA text: the name is the first word after
// '>', the sequence the letters of the following lines (uppercase).
// Not noexcept: the strings allocate.
[[nodiscard]] constexpr auto parse_fasta(std::string_view const text) -> MoleculeResult {
    constexpr auto header_mark = std::string_view{">"};
    auto lines = detail::Lines{text};
    auto line = lines.next();
    while (line.has_value() and detail::is_blank_line(*line)) {
        line = lines.next();
    }
    if (not line.has_value() or not line->starts_with(header_mark)) {
        return detail::failure(MoleculeError::missing_header);
    }
    auto name = std::string{detail::first_word_after(*line, header_mark)};
    auto sequence = std::string{};
    for (line = lines.next(); line.has_value() and not line->starts_with(header_mark); line = lines.next()) {
        detail::append_letters(sequence, *line);
    }
    return detail::success(std::move(name), std::move(sequence));
}

// A file holding only a sequence: its letters (uppercase), named after
// the file. Not noexcept: the strings allocate.
[[nodiscard]] constexpr auto parse_plain_sequence(MoleculeFile const& file) -> MoleculeResult {
    auto sequence = std::string{};
    detail::append_letters(sequence, file.text);
    return detail::success(std::string{file.name}, std::move(sequence));
}

// Any of the three formats, recognised from the first line that is not
// blank (LOCUS: GenBank; '>': FASTA; otherwise a plain sequence). How
// Primer Designer recognises formats was not recorded; this choice
// loads every recorded file as it does. Not noexcept: allocates.
[[nodiscard]] constexpr auto parse_molecule_file(MoleculeFile const& file) -> MoleculeResult {
    auto const text = file.text;
    auto const first = detail::first_content_line(text);
    if (first.has_value() and first->starts_with("LOCUS")) {
        auto genbank = parse_genbank(text);
        if (not genbank.molecule.has_value()) {
            return detail::failure(MoleculeError::missing_origin);
        }
        return detail::success(std::move(genbank.molecule->name), std::move(genbank.molecule->sequence));
    }
    if (first.has_value() and first->starts_with('>')) {
        return parse_fasta(text);
    }
    return parse_plain_sequence(file);
}

}  // namespace primerpicker::io
