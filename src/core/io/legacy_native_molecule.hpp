#pragma once

#include "core/io/genbank.hpp"
#include "core/sequence/topology.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace primerpicker::io {

// The molecule files Primer Designer 2.0 writes (Molecule-Ops ->
// Save/Remove -> S, Create New; tests/legacy/README.md, "Saving and
// creating molecules"): one value a line, DOS line ends:
//   name, size, 0,
//   the number of genes, then per gene: name NUL description, start,
//   end, cw or ccw, Y or N,
//   the number of markers, then per marker: name NUL description,
//   position, R or L,
//   linear molecules only: linear, left end, right end, and two lines
//   ("^" when no end was typed),
//   if there is a sequence: sequence (DNA) or sequenceAA (protein), the
//   whole sequence,
//   if there is a description: description, the text (read before or
//   after the sequence, written after it).
// REGION features are not saved.

enum struct NativeError : std::uint8_t { none, illegal_format };
enum struct GeneDirection : std::uint8_t { clockwise, counterclockwise };  // cw, ccw
enum struct ComplementFlag : std::uint8_t { no, yes };                     // N, Y
enum struct MarkerDirection : std::uint8_t { right, left };                // R, L ("Complement: Y")
enum struct SequenceKind : std::uint8_t { dna, protein };                  // sequence, sequenceAA

struct NativeGene {
    std::string name;
    std::string description;
    std::size_t start;
    std::size_t end;
    GeneDirection direction;  // as read; written from shows_complement
    ComplementFlag flag;
};

struct NativeMarker {
    std::string name;
    std::string description;
    std::size_t position;
    MarkerDirection direction;
};

// the four lines after "linear"
struct NativeEnds {
    std::string left;
    std::string right;
    std::string left_mark;   // "^" when no end was typed
    std::string right_mark;  // the same
};

struct NativeMolecule {
    std::string name;
    std::size_t size = 0;
    std::vector<NativeGene> genes;
    std::vector<NativeMarker> markers;
    std::optional<NativeEnds> linear_ends;  // none: circular
    std::string description;                // empty: not written
    std::optional<std::string> sequence;    // as written (upper case)
    SequenceKind kind = SequenceKind::dna;
};

// C++23 refactoring: replace with std::expected<NativeMolecule, NativeError>
struct NativeResult {
    std::optional<NativeMolecule> molecule;
    NativeError error = NativeError::none;
};

namespace native_detail {

constexpr auto line_end = std::string_view{"\r\n"};
constexpr auto separator = '\0';  // between a feature's name and description
constexpr auto unused_field = std::string_view{"0"};
constexpr auto linear_mark = std::string_view{"linear"};
constexpr auto description_mark = std::string_view{"description"};
constexpr auto dna_mark = std::string_view{"sequence"};
constexpr auto protein_mark = std::string_view{"sequenceAA"};
constexpr auto clockwise_mark = std::string_view{"cw"};
constexpr auto counterclockwise_mark = std::string_view{"ccw"};
constexpr auto yes_mark = std::string_view{"Y"};
constexpr auto no_mark = std::string_view{"N"};
constexpr auto right_mark = std::string_view{"R"};
constexpr auto left_mark = std::string_view{"L"};
constexpr auto no_end_mark = std::string_view{"^"};

constexpr auto failure() -> NativeResult {
    return {.molecule = std::nullopt, .error = NativeError::illegal_format};
}

// a decimal number, digits only (none if it does not fit)
// C++23 refactoring: use std::from_chars (constexpr)
constexpr auto number(std::optional<std::string_view> const line) noexcept -> std::optional<std::size_t> {
    constexpr auto base = std::size_t{10};
    if (not line.has_value() or line->empty()) {
        return std::nullopt;
    }
    auto value = std::size_t{0};
    for (auto const character : *line) {
        if (character < '0' or character > '9') {
            return std::nullopt;
        }
        auto const digit = static_cast<std::size_t>(character - '0');
        if (value > (std::numeric_limits<std::size_t>::max() - digit) / base) {
            return std::nullopt;
        }
        value = (value * base) + digit;
    }
    return value;
}

// "name NUL description"; a line without NUL is shown as both
// (tests/legacy/molecule_save/hand3)
struct Names {
    std::string name;
    std::string description;
};

constexpr auto names(std::string_view const line) -> Names {
    auto const split = line.find(separator);
    if (split == std::string_view::npos) {
        return {.name = std::string{line}, .description = std::string{line}};
    }
    return {.name = std::string{line.substr(0, split)}, .description = std::string{line.substr(split + 1)}};
}

constexpr auto read_gene(detail::Lines& lines) -> std::optional<NativeGene> {
    auto const label = lines.next();
    auto const start = number(lines.next());
    auto const end = number(lines.next());
    auto const direction = lines.next();
    auto const flag = lines.next();
    if (not label.has_value() or not start.has_value() or not end.has_value() or not direction.has_value() or
        not flag.has_value()) {
        return std::nullopt;
    }
    auto [name, description] = names(*label);
    // any other value loads (hand3: R): read as the default
    return NativeGene{.name = std::move(name),
                      .description = std::move(description),
                      .start = *start,
                      .end = *end,
                      .direction = *direction == counterclockwise_mark ? GeneDirection::counterclockwise
                                                                       : GeneDirection::clockwise,
                      .flag = *flag == yes_mark ? ComplementFlag::yes : ComplementFlag::no};
}

constexpr auto read_marker(detail::Lines& lines) -> std::optional<NativeMarker> {
    auto const label = lines.next();
    auto const position = number(lines.next());
    auto const direction = lines.next();
    if (not label.has_value() or not position.has_value() or not direction.has_value()) {
        return std::nullopt;
    }
    auto [name, description] = names(*label);
    return NativeMarker{.name = std::move(name),
                        .description = std::move(description),
                        .position = *position,
                        .direction = *direction == left_mark ? MarkerDirection::left : MarkerDirection::right};
}

constexpr auto read_ends(detail::Lines& lines) -> std::optional<NativeEnds> {
    auto const left = lines.next();
    auto const right = lines.next();
    auto const left_end_mark = lines.next();
    auto const right_end_mark = lines.next();
    if (not left.has_value() or not right.has_value() or not left_end_mark.has_value() or
        not right_end_mark.has_value()) {
        return std::nullopt;
    }
    return NativeEnds{.left = std::string{*left},
                      .right = std::string{*right},
                      .left_mark = std::string{*left_end_mark},
                      .right_mark = std::string{*right_end_mark}};
}

// C++26 refactoring: replace with std::to_string (constexpr)
constexpr auto decimal(std::size_t value) -> std::string {
    constexpr auto base = std::size_t{10};
    auto digits = std::string(1, static_cast<char>('0' + (value % base)));
    value /= base;
    while (value != 0) {
        digits.insert(digits.begin(), static_cast<char>('0' + (value % base)));
        value /= base;
    }
    return digits;
}

constexpr auto is_blank(std::string_view const line) noexcept -> bool {
    return line.find_first_not_of(" \t") == std::string_view::npos;
}

// the sections after the features: the linear ends, then a description
// and a sequence in either order, each optional; then only blank lines
constexpr auto read_tail(detail::Lines& lines, NativeMolecule& molecule) -> bool {
    auto line = lines.next();
    if (line == linear_mark) {
        molecule.linear_ends = read_ends(lines);
        if (not molecule.linear_ends.has_value()) {
            return false;
        }
        line = lines.next();
    }
    auto description_read = false;
    for (; line.has_value(); line = lines.next()) {
        if (line == description_mark and not description_read) {
            auto const text = lines.next();
            if (not text.has_value()) {
                return false;
            }
            molecule.description = std::string{*text};
            description_read = true;
            continue;
        }
        if ((line == dna_mark or line == protein_mark) and not molecule.sequence.has_value()) {
            molecule.kind = line == protein_mark ? SequenceKind::protein : SequenceKind::dna;
            auto const letters = lines.next();
            if (not letters.has_value()) {
                return false;
            }
            molecule.sequence = std::string{*letters};
            continue;
        }
        if (not is_blank(*line)) {
            return false;
        }
    }
    return true;
}

}  // namespace native_detail

// Not noexcept: the strings and vectors allocate.
[[nodiscard]] constexpr auto parse_native_molecule(std::string_view const text) -> NativeResult {
    using namespace native_detail;
    auto lines = detail::Lines{text};
    auto const name = lines.next();
    auto const size = number(lines.next());
    // PRIMER.exe writes 0 here, and rejects other values
    // (tests/legacy/molecule_save/hand1, hand9, handA)
    if (not name.has_value() or not size.has_value() or lines.next() != unused_field) {
        return failure();
    }
    auto molecule = NativeMolecule{.name = std::string{*name},
                                   .size = *size,
                                   .genes = {},
                                   .markers = {},
                                   .linear_ends = std::nullopt,
                                   .description = {},
                                   .sequence = std::nullopt,
                                   .kind = SequenceKind::dna};
    auto const gene_count = number(lines.next());
    if (not gene_count.has_value()) {
        return failure();
    }
    for (auto index = std::size_t{0}; index < *gene_count; ++index) {
        auto gene = read_gene(lines);
        if (not gene.has_value()) {
            return failure();
        }
        molecule.genes.push_back(std::move(*gene));
    }
    auto const marker_count = number(lines.next());
    if (not marker_count.has_value()) {
        return failure();
    }
    for (auto index = std::size_t{0}; index < *marker_count; ++index) {
        auto marker = read_marker(lines);
        if (not marker.has_value()) {
            return failure();
        }
        molecule.markers.push_back(std::move(*marker));
    }
    if (not read_tail(lines, molecule)) {
        return failure();
    }
    return {.molecule = std::move(molecule), .error = NativeError::none};
}

// The features table shows "(C)" for a gene whose flag is Y and whose
// start is before its end, or whose flag is N and whose start is after
// its end (tests/legacy/molecule_save/hand6, hand7, hand8)
[[nodiscard]] constexpr auto shows_complement(NativeGene const& gene) noexcept -> bool {
    return (gene.flag == ComplementFlag::yes) != (gene.start > gene.end);
}

// Not noexcept: the string allocates.
[[nodiscard]] constexpr auto write_native_molecule(NativeMolecule const& molecule) -> std::string {
    using namespace native_detail;
    auto text = std::string{};
    auto const add = [&text](std::string_view const line) {
        text += line;
        text += line_end;
    };
    auto const add_names = [&text, &add](std::string const& name, std::string const& description) {
        text += name;
        text += separator;
        add(description);
    };
    add(molecule.name);
    add(decimal(molecule.size));
    add(unused_field);
    add(decimal(molecule.genes.size()));
    for (auto const& gene : molecule.genes) {
        add_names(gene.name, gene.description);
        add(decimal(gene.start));
        add(decimal(gene.end));
        // ccw when the gene shows (C), whatever was read
        // (tests/legacy/molecule_save/round_zh6)
        add(shows_complement(gene) ? counterclockwise_mark : clockwise_mark);
        add(gene.flag == ComplementFlag::yes ? yes_mark : no_mark);
    }
    add(decimal(molecule.markers.size()));
    for (auto const& marker : molecule.markers) {
        add_names(marker.name, marker.description);
        add(decimal(marker.position));
        add(marker.direction == MarkerDirection::left ? left_mark : right_mark);
    }
    if (molecule.linear_ends.has_value()) {
        add(linear_mark);
        add(molecule.linear_ends->left);
        add(molecule.linear_ends->right);
        add(molecule.linear_ends->left_mark);
        add(molecule.linear_ends->right_mark);
    }
    if (molecule.sequence.has_value()) {
        add(molecule.kind == SequenceKind::protein ? protein_mark : dna_mark);
        add(*molecule.sequence);
    }
    if (not molecule.description.empty()) {
        add(description_mark);
        add(molecule.description);
    }
    return text;
}

// A molecule loaded from a GenBank, FASTA or plain-sequence file, as
// Save writes it: no feature, no description (DEFINITION and FASTA
// descriptions are not kept), no end typed. Not noexcept: copies.
[[nodiscard]] constexpr auto native_from_molecule(Molecule const& molecule) -> NativeMolecule {
    using namespace native_detail;
    auto native = NativeMolecule{.name = molecule.name,
                                 .size = molecule.sequence.size(),
                                 .genes = {},
                                 .markers = {},
                                 .linear_ends = std::nullopt,
                                 .description = {},
                                 .sequence = molecule.sequence,
                                 .kind = SequenceKind::dna};
    if (molecule.topology == sequence::Topology::linear) {
        native.linear_ends = NativeEnds{
            .left = {}, .right = {}, .left_mark = std::string{no_end_mark}, .right_mark = std::string{no_end_mark}};
    }
    return native;
}


}  // namespace primerpicker::io
