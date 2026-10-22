#include "core/report/legacy_report.hpp"

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_duplex.hpp"
#include "core/analysis/legacy_homology.hpp"
#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_pair.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/search/legacy_match.hpp"
#include "core/sequence/genetic_code.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <optional>
#include <cstddef>
#include <format>
#include <iterator>
#include <span>
#include <string>
#include <string_view>
#include <utility>

namespace primerpicker::report {

namespace {

using analysis::PrimerPair;
using sequence::Nucleotide;

constexpr auto title_width = std::size_t{54};  // "Report Date:" starts in this column
constexpr auto lines_after_title = std::size_t{5};
constexpr auto cp437_degree_sign = std::string_view{"\xF8"};  // the reports use code page 437
constexpr auto molecule_name_width = std::size_t{34};
constexpr auto molecule_size_width = std::size_t{5};  // right-aligned, ending in column 56
constexpr auto years_per_century = 100;
constexpr auto criteria_name_width = std::size_t{17};
constexpr auto criteria_actual_end = std::size_t{56};  // right edge of the "Actual" column
constexpr auto verdict_gap = std::string_view{"           "};
constexpr auto duplex_column = std::ptrdiff_t{13};  // the first primer's 5' base
constexpr auto duplex_indent = std::string_view{"          "};
constexpr auto five_prime_label = std::string_view{"5'-"};
constexpr auto three_prime_label = std::string_view{"3'-"};
constexpr auto three_prime_end = std::string_view{"-3'"};
constexpr auto five_prime_end = std::string_view{"-5'"};
constexpr auto page_break = std::string_view{"\r\f"};  // before page 2
constexpr auto rule = std::string_view{"---------------------------------------------------------------------------"};

auto line(std::string_view const text) -> std::string {
    return std::string{text} + std::string{legacy_line_end};
}

auto blank_lines(std::size_t const count) -> std::string {
    auto text = std::string{};
    for (auto index = std::size_t{0}; index < count; ++index) {
        text += legacy_line_end;
    }
    return text;
}

// The width that right-aligns TEXT_SIZE characters to end before column
// END, USED columns being already written; at least one space when they
// do not fit (long settings or sites, never recorded from PRIMER.exe)
constexpr auto right_aligned_width(std::size_t const used, std::size_t const text_size,
                                   std::size_t const end) noexcept -> std::size_t {
    return used + text_size < end ? end - used : text_size + 1;
}

auto letters(std::span<Nucleotide const> const primer) -> std::string {
    auto text = std::string{};
    std::ranges::transform(primer, std::back_inserter(text), sequence::to_letter);
    return text;
}

// "EXTENSION-PRIMER" when there is a 5' extension, as PRIMER.exe's
// reports write it
auto letters(analysis::LegacyPrimer const& primer) -> std::string {
    auto const sequence = primer.sequence;
    if (primer.extension_length == 0) {
        return letters(sequence);
    }
    return letters(sequence.first(primer.extension_length)) + "-" + letters(sequence.subspan(primer.extension_length));
}

auto title_line(std::string_view const title, std::chrono::year_month_day const& date) -> std::string {
    return line(std::format("{:<{}}Report Date: {:02}/{:02}/{:02}", title, title_width,
                            static_cast<unsigned>(date.month()), static_cast<unsigned>(date.day()),
                            static_cast<int>(date.year()) % years_per_century));
}

auto verdict_word(bool const meets) -> std::string_view {
    return meets ? "YES" : "NO";
}

// One row of the criteria table
struct Criterion {
    std::string_view name;
    std::string setting;
    std::string actual;
    bool meets;
};

// The setting after the criterion's name, the actual value right-aligned
// on criteria_actual_end, then the verdict
auto criterion_line(Criterion const& criterion) -> std::string {
    auto const left = std::format("{:<{}}{}", criterion.name, criteria_name_width, criterion.setting);
    return line(std::format("{}{:>{}}{}{}", left, criterion.actual,
                            right_aligned_width(left.size(), criterion.actual.size(), criteria_actual_end), verdict_gap,
                            verdict_word(criterion.meets)));
}

auto range_setting(int const minimum, int const maximum) -> std::string {
    return std::format("Min {:4},    Max {:4}", minimum, maximum);
}

auto criteria_rows(analysis::LegacyActualValues const& actual, analysis::LegacyCriteria const& criteria,
                   analysis::LegacyVerdict const& verdict) -> std::string {
    auto rows = criterion_line({.name = "% GC",
                                .setting = range_setting(criteria.gc_min_percent, criteria.gc_max_percent),
                                .actual = std::to_string(actual.gc_percent),
                                .meets = verdict.gc});
    rows += legacy_line_end;
    rows += criterion_line({.name = "Tm C",
                            .setting = range_setting(criteria.tm_min_celsius, criteria.tm_max_celsius),
                            .actual = std::to_string(actual.tm_celsius),
                            .meets = verdict.tm});
    rows += legacy_line_end;
    // as on the screen, the hairpin energy is not in the "Actual" column
    rows += criterion_line({.name = "No Hairpins",
                            .setting = std::format("Energy cutoff {:4.1f} kcal", criteria.hairpin_cutoff.kcal_per_mol),
                            .actual = "--",
                            .meets = verdict.no_hairpin});
    rows += legacy_line_end;
    rows += criterion_line(
        {.name = "No 3' Dimers",
         .setting = std::format("Reject >= {:2} matches at 3' end", criteria.reject_dimer_3prime_matches),
         .actual = std::to_string(actual.dimer_3prime_matches),
         .meets = verdict.no_3prime_dimer});
    rows += criterion_line(
        {.name = "No Dimers",
         .setting = std::format("Reject >= {:2} adj homol bases", criteria.reject_dimer_adjacent_bases),
         .actual = std::to_string(actual.dimer_adjacent_bases),
         .meets = verdict.no_dimer});
    rows += legacy_line_end;
    rows += criterion_line({.name = "No Runs",
                            .setting = std::format("Reject >= {:2} base runs", criteria.reject_run_length),
                            .actual = std::to_string(actual.longest_run),
                            .meets = verdict.no_runs});
    rows += criterion_line({.name = "No 3'GC runs",
                            .setting = std::format("Reject >= {:2} G or C at 3' end", criteria.reject_gc_run_3prime),
                            .actual = std::to_string(actual.gc_run_3prime),
                            .meets = verdict.no_3prime_gc_run});
    return rows;
}

// The second primer 3' to 5', its 3' base in COLUMN: the "3'-" label is
// left out when it does not fit left of that base, and bases left of
// column 0 are cut
auto second_primer_line(std::span<Nucleotide const> const second, std::ptrdiff_t const column) -> std::string {
    auto reversed = letters(second);
    std::ranges::reverse(reversed);
    auto const label = static_cast<std::ptrdiff_t>(three_prime_label.size());
    if (column >= label) {
        return std::string(static_cast<std::size_t>(column - label), ' ') + std::string{three_prime_label} +
               reversed + std::string{five_prime_end};
    }
    if (column >= 0) {
        return std::string(static_cast<std::size_t>(column), ' ') + reversed + std::string{five_prime_end};
    }
    auto const hidden = static_cast<std::size_t>(-column);
    return (hidden < reversed.size() ? reversed.substr(hidden) : std::string{}) + std::string{five_prime_end};
}

// A duplex without pairs: PRIMER.exe writes lines of tens of thousands
// of spaces there (tests/legacy/README.md, "Duplexes without pairs").
// Decided with the developer (2026-10-20): the partner goes right under
// the primer, its "3'-" under the "5'-", below an empty marks line.
constexpr auto no_pairs_3prime_column = static_cast<std::ptrdiff_t>(duplex_indent.size() + three_prime_label.size());

auto duplex_lines(PrimerPair const& pair_of_primers, analysis::DuplexPlacement const& placement) -> std::string {
    auto text = line(std::string{duplex_indent} + std::string{five_prime_label} + letters(pair_of_primers.first) +
                     std::string{three_prime_end});
    if (placement.count == 0) {
        text += legacy_line_end;
        text += line(second_primer_line(pair_of_primers.second, no_pairs_3prime_column));
        return text;
    }
    auto const marks_column = duplex_column + placement.marks_offset;
    assert(marks_column >= 0);
    text += line(std::string(static_cast<std::size_t>(marks_column), ' ') + std::string(placement.count, ':'));
    text += line(second_primer_line(pair_of_primers.second, duplex_column + placement.second_3prime_offset));
    return text;
}


// Primer Pair report: values right-aligned on these columns (inclusive),
// each followed by its flag ('!' or a space)
constexpr auto pair_label_column = std::size_t{14};
constexpr auto pair_first_value_end = std::size_t{40};  // one past the value of A
constexpr auto pair_value_step = std::size_t{9};       // to the next primer
constexpr auto balance_comment_gap = std::string_view{"      "};  // after the value of B
constexpr auto caution_comment_gap = std::string_view{"     "};
constexpr auto pair_dimer_note_gap = std::string_view{"   "};

auto flag(bool const flagged) -> char {
    return flagged ? '!' : ' ';
}

// LABEL, then each value with its flag, as the rows of the Primer Pair
// report
struct FlaggedValue {
    std::string text;
    bool flagged;
};

auto pair_row(std::string_view const label, std::span<FlaggedValue const> const values) -> std::string {
    auto row = std::format("{:{}}{}", "", pair_label_column, label);
    auto end = pair_first_value_end;
    for (auto const& value : values) {
        row += std::format("{:>{}}{}", value.text, right_aligned_width(row.size(), value.text.size(), end),
                           flag(value.flagged));
        end += pair_value_step;
    }
    return row;
}

// The criteria in short form (Primer Pair and Find Pair reports)
auto short_criteria(analysis::LegacyCriteria const& criteria) -> std::string {
    auto text = line(std::format("Criteria:     % GC:  Min {:3}, Max {:3};    Tm C:  Min {:3}, Max {:3}",
                                 criteria.gc_min_percent, criteria.gc_max_percent, criteria.tm_min_celsius,
                                 criteria.tm_max_celsius));
    text += line(std::format("              Hairpin cutoff:{:4.1f} kcal;", criteria.hairpin_cutoff.kcal_per_mol));
    text += line(std::format("              Reject >= {} matches at 3' end;  >= {} adj homol bases",
                             criteria.reject_dimer_3prime_matches, criteria.reject_dimer_adjacent_bases));
    text += line(std::format("              Reject >= {} base runs;  >= {} G/C at 3' end", criteria.reject_run_length,
                             criteria.reject_gc_run_3prime));
    return text;
}

// A false-priming Tm of the Find Pair list: '-' when there is no site
auto site_tm(std::optional<int> const tm_celsius, std::size_t const width) -> std::string {
    return tm_celsius.has_value() ? std::format("{:>{}}", *tm_celsius, width) : std::string{"-"};
}

auto find_pair_row(search::PairCandidate const& pair) -> std::string {
    constexpr auto tm_b_width = std::size_t{2};
    return line(std::format("{:>8}{:>10}{:>8} / {:<2}{:>6} / {:<2}{:>8}{:>5}{:>9} / {}", pair.position_a,
                            pair.position_b, pair.primer_a.gc_percent, pair.primer_b.gc_percent,
                            pair.primer_a.tm_celsius, pair.primer_b.tm_celsius, pair.max_dimer_3prime_matches,
                            pair.max_dimer_adjacent_bases, site_tm(pair.false_priming_tm_a, 0),
                            site_tm(pair.false_priming_tm_b, tm_b_width)));
}

auto primer_header(char const name) -> std::string {
    if (name == 'A') {
        return "-----  P R I M E R   A  ---------------------------------------------------";
    }
    return "-----------------------------------------------  P R I M E R   B  ---------";
}

auto primer_summary(char const name, analysis::LegacyPrimer const& extended, analysis::LegacyCriteria const& criteria)
    -> std::string {
    auto const primer = extended.sequence;
    auto const actual = analysis::legacy_actual_values(extended);
    auto const verdict = analysis::evaluate(actual, criteria);
    auto const count = [primer](Nucleotide const nucleotide) { return std::ranges::count(primer, nucleotide); };
    auto text = line(primer_header(name));
    text += legacy_line_end;
    text += line(std::format("Sequence:     5'-{}-3'", letters(extended)));
    text += legacy_line_end;
    text += line(std::format("              Length: {:3};", primer.size()));
    text += line(std::format("              % GC:{:3};     Tm C:{:3};    Composition:  A:{}  C:{}  G:{}  T:{}",
                             actual.gc_percent, actual.tm_celsius, count(Nucleotide::adenine),
                             count(Nucleotide::cytosine), count(Nucleotide::guanine), count(Nucleotide::thymine)));
    text += legacy_line_end;
    text += short_criteria(criteria);
    text += line(std::format("              Meets criteria selected:  {}", verdict_word(verdict.all())));
    text += legacy_line_end;
    return text;
}

// A row of the primer balance: both values, then "N base difference",
// "N % difference" or "N C difference"
struct Balance {
    std::string_view label;
    int value_a;
    int value_b;
    int difference;
    std::string_view unit;
};

// (no comment when the difference is zero)
auto balance_row(Balance const& balance) -> std::string {
    auto const row = std::format("{:{}}{}", "", pair_label_column, balance.label);
    auto text = std::format("{}{:>{}}{:>{}}", row, balance.value_a, pair_first_value_end - row.size(), balance.value_b,
                            pair_value_step);
    if (balance.difference != 0) {
        text += std::format("{}{}{} difference", balance_comment_gap, balance.difference, balance.unit);
    }
    return line(text);
}

auto caution_row(std::string_view const label, FlaggedValue const& first, FlaggedValue const& second) -> std::string {
    return line(pair_row(label, std::array{first, second}));
}

// "OK" or "!" for each primer, without flags, then "Hairpin(s) found"
// when either primer has one
auto hairpin_row(bool const hairpin_a, bool const hairpin_b) -> std::string {
    auto const text = [](bool const hairpin) { return hairpin ? "!" : "OK"; };
    auto const row = std::format("{:{}}{}", "", pair_label_column, "Hairpins");
    auto hairpins = std::format("{}{:>{}}{:>{}}", row, text(hairpin_a), pair_first_value_end - row.size(),
                                text(hairpin_b), pair_value_step);
    if (hairpin_a or hairpin_b) {
        hairpins += std::format("{}Hairpin(s) found", caution_comment_gap);
    }
    return line(hairpins);
}

// Multi Primer report: text placed at columns of a line, which then
// loses its trailing spaces
class Canvas {
public:
    auto put(std::size_t const column, std::string_view const text) -> Canvas& {
        if (text_.size() < column + text.size()) {
            text_.resize(column + text.size(), ' ');
        }
        text_.replace(column, text.size(), text);
        return *this;
    }

    // TEXT ending in column LAST
    auto put_right(std::size_t const last, std::string_view const text) -> Canvas& {
        assert(text.size() <= last + 1);
        return put(last + 1 - text.size(), text);
    }

    [[nodiscard]] auto line() const -> std::string {
        auto const end = text_.find_last_not_of(' ');
        return report::line(end == std::string::npos ? std::string_view{} : std::string_view{text_}.substr(0, end + 1));
    }

private:
    std::string text_;
};

constexpr auto multi_names = std::string_view{"ABCD"};
constexpr auto multi_label_column = std::size_t{6};
constexpr auto multi_first_value_last = std::size_t{48};  // last column of primer A's values
constexpr auto multi_value_step = std::size_t{8};
constexpr auto matrix_cell_step = std::size_t{6};
constexpr auto matrix_row_name_column = std::size_t{10};
constexpr auto matrix_first_value_last = std::size_t{14};  // "3' end" matrix, column A
constexpr auto matrix_dots_column = std::size_t{13};
constexpr auto matrix_dots_shorter = std::size_t{2};       // than the cells they cover
constexpr auto matrix_offset = std::size_t{37};            // from the "3' end" matrix to the "Any" one

// A field's section: its primer, or "none" for an empty field
auto multi_primer_section(std::size_t const index, analysis::LegacyMultiPrimerFields const& fields) -> std::string {
    auto text = line(std::format("-----  P R I M E R   {}  ---------------------------------------------------",
                                 multi_names.at(index)));
    text += legacy_line_end;
    auto const& field = fields.at(index);
    if (not field.has_value()) {
        return text + line("     none") + std::string{legacy_line_end};
    }
    auto const& primer = *field;
    text += line(std::format("Sequence:     5'-{}-3'", letters(primer)));
    text += legacy_line_end;
    text += line("              Source molecule:  none");
    text += line(std::format("              Length: {:3};", primer.sequence.size()));
    text += legacy_line_end;
    return text;
}

// A row of "Primer Balance and Cautions": one value per primer, no
// flag; an empty field leaves its column blank (trailing spaces kept)
template <typename Value>
auto multi_balance_row(std::string_view const label, analysis::LegacyMultiPrimer const& analysis, Value const value)
    -> std::string {
    auto row = std::format("{:{}}{:<{}}", "", multi_label_column, label,
                           multi_first_value_last + 1 - multi_value_step - multi_label_column);
    for (auto const& column : analysis.columns) {
        row += std::format("{:>{}}", column.has_value() ? value(*column) : std::string{}, multi_value_step);
    }
    return line(row);
}

// Row ROW of both dimer matrices: dots below the diagonal, values with
// their '!' flag on and above it; an empty field's row has neither
auto matrix_row(std::size_t const row, analysis::LegacyMultiPrimer const& analysis) -> std::string {
    constexpr auto three_prime_labels = std::array{std::string_view{"3'"}, std::string_view{"end"}};
    constexpr auto label_column = std::size_t{3};
    constexpr auto any_label_column = std::size_t{40};
    auto canvas = Canvas{};
    if (row < three_prime_labels.size()) {
        canvas.put(label_column, three_prime_labels.at(row));
    }
    if (row == 0) {
        canvas.put(any_label_column, "Any");
    }
    auto const name = std::string(1, multi_names.at(row));
    for (auto const offset : {std::size_t{0}, matrix_offset}) {
        canvas.put(matrix_row_name_column + offset, name);
        if (row > 0 and analysis.columns.at(row).has_value()) {
            canvas.put(matrix_dots_column + offset, std::string((row * matrix_cell_step) - matrix_dots_shorter, '.'));
        }
    }
    for (auto column = row; column < analysis::legacy_multi_primer_fields; ++column) {
        auto const& cell = analysis.dimers.at(row).at(column);
        if (not cell.has_value()) {
            continue;
        }
        auto const last = matrix_first_value_last + (column * matrix_cell_step);
        auto const flagged = [](std::size_t const value, bool const is_flagged) {
            return std::to_string(value) + (is_flagged ? "!" : "");
        };
        canvas.put(last + 1 - std::to_string(cell->dimer_3prime_matches).size(),
                   flagged(cell->dimer_3prime_matches, cell->flags.three_prime));
        canvas.put(last + matrix_offset + 1 - std::to_string(cell->dimer_adjacent_bases).size(),
                   flagged(cell->dimer_adjacent_bases, cell->flags.adjacent));
    }
    return canvas.line();
}

// A duplex of page 2, under its name ("A:B")
auto named_duplex(std::string_view const name, PrimerPair const& pair_of_primers,
                  analysis::DuplexPlacement const& placement) -> std::string {
    return line(name) + duplex_lines(pair_of_primers, placement) + blank_lines(2);
}

// The text map of the Uniqueness and False Priming reports: the molecule
// drawn with one character per molecule_size / 73 + 1 bases
constexpr auto map_scale_divisor = std::size_t{73};
constexpr auto map_site_margin = std::size_t{3};  // the column of position 1
constexpr auto homology_bases_width = std::size_t{19};  // dots, then the matching bases
constexpr auto homology_position_end = std::size_t{53};  // one past the "Bp Pos" column
constexpr auto homology_bases_end = std::size_t{34};      // one past the bases' last column
constexpr auto long_site_prefix = std::string_view{" 5'"};
constexpr auto homology_tm_width = std::size_t{10};

auto bases_per_character(std::size_t const molecule_size) -> std::size_t {
    return (molecule_size / map_scale_divisor) + 1;
}

// A line of site marks ('>' top strand, '<' complement), or a single
// space when there is none
auto map_marks(std::span<analysis::HomologySite const> const sites, sequence::Strand const strand,
               std::size_t const scale) -> std::string {
    auto canvas = Canvas{};
    auto const mark = strand == sequence::Strand::top ? std::string_view{">"} : std::string_view{"<"};
    for (auto const& site : sites) {
        if (site.strand == strand) {
            canvas.put(map_site_margin + ((site.position - 1) / scale), mark);
        }
    }
    auto marks = canvas.line();
    return marks == line("") ? line(" ") : marks;
}

auto homology_map(HomologyReport const& homology) -> std::string {
    auto const scale = bases_per_character(homology.molecule_size);
    auto const primer_column = homology.position / scale;
    auto text = line(std::string(primer_column > 0 ? primer_column - 1 : 0, ' ') + "P>");
    text += line(std::string(homology.molecule_size / scale, '='));
    text += map_marks(homology.sites, sequence::Strand::top, scale);
    text += map_marks(homology.sites, sequence::Strand::complement, scale);
    return text;
}

// One site: the matching bases right-aligned after dots, "....3'" for a
// uniqueness site (the match stops before the 3' end), the position,
// "(C)" on the complement strand, the Tm
auto homology_row(HomologyReport const& homology, analysis::HomologySite const& site) -> std::string {
    auto const bases = letters(homology.primer.subspan(site.primer_start, site.length));
    auto const three_prime = homology.kind == analysis::HomologyKind::uniqueness ? std::string_view{"....3'"}
                                                                              : std::string_view{"      "};
    // a site longer than the field: " 5'" in column 1, the bases still
    // ending in the field's last column (20 and 21 bases recorded)
    auto row = bases.size() <= homology_bases_width
                   ? std::format("             5'{:.>{}}{}", bases, homology_bases_width, three_prime)
                   : std::format(" 5'{:>{}}{}", bases, right_aligned_width(long_site_prefix.size(), bases.size(),
                                                                         homology_bases_end),
                                 three_prime);
    auto const position = std::to_string(site.position);
    row += std::format("{:>{}}", position, right_aligned_width(row.size(), position.size(), homology_position_end));
    row += site.strand == sequence::Strand::complement ? " (C)" : "    ";
    return line(std::format("{}{:>{}}   ", row, site.tm_celsius, homology_tm_width));
}

// Match Active Primer report rows: values right-aligned (primer A's
// position ends in column 5, primer B's in column 7)
auto match_values(analysis::LegacyActualValues const& actual) -> std::string {
    return std::format("{:>8}{:>6}{:>8}{:>5}", actual.gc_percent, actual.tm_celsius, actual.dimer_3prime_matches,
                       actual.dimer_adjacent_bases);
}

auto match_runs(analysis::LegacyActualValues const& actual) -> std::string {
    return std::format("{:>8}{:>6}", actual.longest_run, actual.gc_run_3prime);
}

// Whether primer A's row marks the values that fail their criterion
enum struct PrimerAFlags : std::uint8_t { shown, hidden };

// Primer A's row: a value that fails its criterion is followed by '!'
// (in the column after it)
auto match_primer_a_row(std::size_t const position, analysis::LegacyActualValues const& actual,
                        analysis::LegacyCriteria const& criteria, PrimerAFlags const shown_flags) -> std::string {
    constexpr auto no_distance = std::string_view{"            "};
    constexpr auto no_inter_dimers = std::string_view{"             "};
    // the last column of each value
    constexpr auto gc_end = std::size_t{25};
    constexpr auto tm_end = std::size_t{31};
    constexpr auto dimer_3prime_end = std::size_t{39};
    constexpr auto dimer_any_end = std::size_t{44};
    constexpr auto runs_end = std::size_t{65};
    constexpr auto gc_run_end = std::size_t{71};
    auto row = std::format("{:>6}{}{}{}{}", position, no_distance, match_values(actual), no_inter_dimers,
                           match_runs(actual));
    auto const verdict = analysis::evaluate(actual, criteria);
    auto const flags = std::array{std::pair{gc_end, verdict.gc}, std::pair{tm_end, verdict.tm},
                                  std::pair{dimer_3prime_end, verdict.no_3prime_dimer},
                                  std::pair{dimer_any_end, verdict.no_dimer}, std::pair{runs_end, verdict.no_runs},
                                  std::pair{gc_run_end, verdict.no_3prime_gc_run}};
    for (auto const& [end, meets] : flags) {
        if (not meets and shown_flags == PrimerAFlags::shown) {
            row.resize(std::max(row.size(), end + 2), ' ');
            row[end + 1] = '!';
        }
    }
    return line(row);
}

// The offsets as searched (complement strand): the start cut so that the
// first primer fits (printed as an unsigned 16-bit number when the
// target is shorter than a primer: 65526 for target 10), the end cut at
// the molecule's end
auto match_offsets(MatchReport const& match) -> std::pair<std::size_t, std::size_t> {
    constexpr auto sixteen_bits = std::size_t{65536};
    auto const& settings = match.settings;
    auto const start = settings.target >= settings.length + settings.start_offset
                           ? settings.start_offset
                           : (settings.target + sixteen_bits - settings.length) % sixteen_bits;
    auto const end = std::min(settings.end_offset,
                              match.molecule_size > settings.target ? match.molecule_size - settings.target : 0);
    return {start, end};
}

auto match_row(search::MatchCandidate const& candidate) -> std::string {
    return line(std::format("{:>8}{:>10}{}{:>8}{:>5}{}", candidate.position, candidate.distance,
                            match_values(candidate.actual), candidate.inter_dimer_3prime_matches,
                            candidate.inter_dimer_adjacent_bases, match_runs(candidate.actual)));
}

}  // namespace

namespace {

// What the One Primer report of a mutated primer adds before the
// criteria: the source molecule, "Mutated", the position
struct MutatedHeader {
    std::string_view molecule_name;
    std::size_t molecule_size;
    std::size_t position;
};

auto inspect_report(analysis::LegacyPrimer const& primer, analysis::LegacyCriteria const& criteria,
                    std::chrono::year_month_day const& date, std::optional<MutatedHeader> const& mutated)
    -> std::string {
    auto const actual = analysis::legacy_actual_values(primer);
    auto const verdict = analysis::evaluate(actual, criteria);
    auto const& sequence = primer.sequence;
    auto const extension = sequence.first(primer.extension_length);
    auto const self = PrimerPair{.first = sequence, .second = sequence};
    auto const count = [sequence](Nucleotide const nucleotide) { return std::ranges::count(sequence, nucleotide); };

    auto report = title_line("Primer Design Report:  One Primer", date);
    if (mutated.has_value()) {
        report += blank_lines(2);
        report += line(std::format("Source molecule:  {:<{}}{:>{}} bps", mutated->molecule_name, molecule_name_width,
                                   mutated->molecule_size, molecule_size_width));
        report += blank_lines(3);
    } else {
        report += blank_lines(lines_after_title);
    }
    report += line("-----  P R I M E R   I N F O  --------------------------------------------");
    if (mutated.has_value()) {
        report += legacy_line_end;
        report += line("Mutated");
    } else {
        report += blank_lines(2);
    }
    // the extension, then a dash, then the primer
    auto const shown = extension.empty() ? letters(sequence)
                                         : letters(extension) + "-" + letters(sequence.subspan(extension.size()));
    report += line(std::format("Sequence:     5'-{}-3'", shown));
    report += legacy_line_end;
    if (mutated.has_value()) {
        constexpr auto position_width = std::size_t{5};
        report += line(std::format("              Length: {:3};    Bp pos: {:{}}", sequence.size(), mutated->position,
                                   position_width));
    } else {
        report += line(std::format("              Length: {:3};", sequence.size()));
    }
    report += line(std::format("              5' Extension:   {}",
                               extension.empty() ? std::string{"--none--"} : letters(extension)));
    report += line(std::format("              Composition:    A:{:2}  C:{:2}  G:{:2}  T:{:2}",
                               count(Nucleotide::adenine), count(Nucleotide::cytosine), count(Nucleotide::guanine),
                               count(Nucleotide::thymine)));
    report += legacy_line_end;
    report += line(rule);
    report += blank_lines(3);
    report += line("PRIMER EVALUATION AND CRITERIA SETTINGS:");
    report += blank_lines(2);
    report += line("Criteria         Criteria Setting                  Actual    Meets Criteria");
    report += line("-------------    ----------------------------      ------    --------------");
    report += legacy_line_end;
    report += criteria_rows(actual, criteria, verdict);
    report += blank_lines(4);
    report += line("LARGEST DUPLEXES REPORTED IN PRIMER DIMERS (above)");
    report += blank_lines(2);
    report += line("     Primer-dimers with 3'end involvement:");
    report += legacy_line_end;
    report += duplex_lines(self, analysis::legacy_3prime_duplex(self));
    report += blank_lines(2);
    report += line("     Other Primer-Dimers:");
    report += legacy_line_end;
    report += duplex_lines(self, analysis::legacy_longest_duplex(self));
    return report;
}

constexpr auto mutation_column = std::string_view{"                                   "};  // the sequences' column

// a translation row: each amino acid under the third base of its codon,
// for the codons that start at base FIRST (0-based) and every third
// base after it, inside the primer; as long as the primer, or ending
// with the last amino acid
auto translation_row(std::span<Nucleotide const> const primer, std::size_t const first) -> std::string {
    constexpr auto codon = std::size_t{3};
    auto row = std::string(primer.size(), ' ');
    auto last = std::size_t{0};
    for (auto start = first; start + codon <= primer.size(); start += codon) {
        row.at(start + 2) = sequence::translate(primer[start], primer[start + 1], primer[start + 2]);
        last = start + 2;
    }
    // the last frame's row stops at its last amino acid
    if (first == 2) {
        row.resize(last + 1);
    }
    return row;
}

auto translations(std::string_view const label, std::span<Nucleotide const> const primer) -> std::string {
    auto text = line(std::format("          {} sequence:       {}", label, letters(primer)));
    for (auto const first : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
        text += line(std::string{mutation_column} + translation_row(primer, first));
    }
    return text;
}

}  // namespace

auto legacy_inspect_report(analysis::LegacyPrimer const& primer, analysis::LegacyCriteria const& criteria,
                           std::chrono::year_month_day const& date) -> std::string {
    return inspect_report(primer, criteria, date, std::nullopt);
}

auto legacy_mutated_inspect_report(MutatedPrimerReport const& mutated, analysis::LegacyCriteria const& criteria,
                                   std::chrono::year_month_day const& date) -> std::string {
    assert(mutated.original.size() == mutated.mutated.size());
    auto report = inspect_report({.sequence = mutated.mutated, .extension_length = 0}, criteria, date,
                                 MutatedHeader{.molecule_name = mutated.molecule_name,
                                               .molecule_size = mutated.molecule_size,
                                               .position = mutated.position});
    report += page_break;
    report += line("PRIMER MUTAGENESIS INFORMATION");
    report += blank_lines(2);
    report += line("     Modifications made:");
    report += blank_lines(2);
    report += line(std::format("          Original sequence:       {}", letters(mutated.original)));
    auto marks = std::string(mutated.original.size(), ' ');
    auto changed = std::string(mutated.original.size(), ' ');
    for (auto index = std::size_t{0}; index < mutated.original.size(); ++index) {
        if (mutated.original[index] != mutated.mutated[index]) {
            marks.at(index) = '|';
            changed.at(index) = sequence::to_letter(mutated.mutated[index]);
        }
    }
    report += line(std::string{mutation_column} + marks);
    report += line("            Base(s) changed:       " + changed);
    report += blank_lines(3);
    report += line("     Translation:");
    report += blank_lines(2);
    report += translations("Original", mutated.original);
    report += blank_lines(3);
    report += translations("Modified", mutated.mutated);
    return report;
}


auto legacy_pair_report(PrimerPair const& primers, analysis::LegacyCriteria const& criteria,
                        std::chrono::year_month_day const& date) -> std::string {
    return legacy_pair_report(analysis::LegacyPrimerPair{.first = {.sequence = primers.first, .extension_length = 0},
                                                         .second = {.sequence = primers.second, .extension_length = 0}},
                              criteria, date);
}

auto legacy_pair_report(analysis::LegacyPrimerPair const& extended, analysis::LegacyCriteria const& criteria,
                        std::chrono::year_month_day const& date) -> std::string {
    auto const primers = PrimerPair{.first = extended.first.sequence, .second = extended.second.sequence};
    auto const& [primer_a, primer_b] = primers;
    auto const actual_a = analysis::legacy_actual_values(extended.first);
    auto const actual_b = analysis::legacy_actual_values(extended.second);
    auto const cautions_a = analysis::legacy_pair_cautions(actual_a, criteria);
    auto const cautions_b = analysis::legacy_pair_cautions(actual_b, criteria);
    auto const values_ab = analysis::legacy_pair_values(extended);
    auto const flags = [&criteria](std::size_t const three_prime, std::size_t const adjacent) {
        return analysis::legacy_dimer_flags(three_prime, adjacent, criteria);
    };
    auto const flags_aa = flags(actual_a.dimer_3prime_matches, actual_a.dimer_adjacent_bases);
    auto const flags_bb = flags(actual_b.dimer_3prime_matches, actual_b.dimer_adjacent_bases);
    auto const flags_ab = flags(values_ab.dimer_3prime_matches, values_ab.dimer_adjacent_bases);
    auto const self_a = PrimerPair{.first = primer_a, .second = primer_a};
    auto const self_b = PrimerPair{.first = primer_b, .second = primer_b};
    auto const counted = [](std::size_t const value, bool const flagged) {
        return FlaggedValue{.text = std::to_string(value), .flagged = flagged};
    };

    auto report = title_line("Primer Design Report:  Primer Pair", date);
    report += blank_lines(lines_after_title + 1);
    report += primer_summary('A', extended.first, criteria);
    report += primer_summary('B', extended.second, criteria);
    report += line(rule);
    report += legacy_line_end;
    report += line("Primer Balance:                      -A-      -B-      Comment");
    report += legacy_line_end;
    report += balance_row({.label = "Length",
                           .value_a = static_cast<int>(primer_a.size()),
                           .value_b = static_cast<int>(primer_b.size()),
                           .difference = values_ab.length_difference,
                           .unit = " base"});
    report += balance_row({.label = "% GC",
                           .value_a = actual_a.gc_percent,
                           .value_b = actual_b.gc_percent,
                           .difference = values_ab.gc_difference_percent,
                           .unit = " %"});
    report += balance_row({.label = "Tm C",
                           .value_a = actual_a.tm_celsius,
                           .value_b = actual_b.tm_celsius,
                           .difference = values_ab.tm_difference_celsius,
                           .unit = " C"});
    report += legacy_line_end;
    report += line("Primer Cautions:");
    report += legacy_line_end;
    report += caution_row("Runs of G/C at 3' end", counted(actual_a.gc_run_3prime, cautions_a.gc_run_3prime),
                          counted(actual_b.gc_run_3prime, cautions_b.gc_run_3prime));
    report += caution_row("Runs of bases (any)", counted(actual_a.longest_run, cautions_a.runs),
                          counted(actual_b.longest_run, cautions_b.runs));
    report += hairpin_row(cautions_a.hairpin, cautions_b.hairpin);
    report += legacy_line_end;
    report += line("Primer Dimers:");
    report += line("                                      A:A      B:B      A:B");
    report += line(pair_row("3'end matches", std::array{counted(actual_a.dimer_3prime_matches, flags_aa.three_prime),
                                                         counted(actual_b.dimer_3prime_matches, flags_bb.three_prime),
                                                         counted(values_ab.dimer_3prime_matches, flags_ab.three_prime)}) +
                   std::string{pair_dimer_note_gap} + "(see page 2");
    report += line(pair_row("Any adj homol bases",
                            std::array{counted(actual_a.dimer_adjacent_bases, flags_aa.adjacent),
                                       counted(actual_b.dimer_adjacent_bases, flags_bb.adjacent),
                                       counted(values_ab.dimer_adjacent_bases, flags_ab.adjacent)}) +
                   std::string{pair_dimer_note_gap} + "for duplexes)");
    report += legacy_line_end;
    report += page_break;
    report += title_line("Primer Design Report:  Primer Pair  -- Page 2", date);
    report += blank_lines(2);
    report += line("Largest Duplexes reported in Primer Interactions");
    report += legacy_line_end;
    report += line(rule);
    report += legacy_line_end;
    report += line("Primer-Dimers with 3'end involvement:");
    report += blank_lines(3);
    report += named_duplex("A:A", self_a, analysis::legacy_3prime_duplex(self_a));
    report += named_duplex("B:B", self_b, analysis::legacy_3prime_duplex(self_b));
    report += named_duplex("A:B", primers, analysis::legacy_3prime_duplex(primers));
    report += line(rule);
    report += legacy_line_end;
    report += line("Other Primer-Dimers:");
    report += blank_lines(3);
    report += named_duplex("A:A", self_a, analysis::legacy_longest_duplex(self_a));
    report += named_duplex("B:B", self_b, analysis::legacy_longest_duplex(self_b));
    report += named_duplex("A:B", primers, analysis::legacy_longest_duplex(primers));
    report += line(rule);
    return report;
}


auto legacy_multi_primer_report(analysis::MultiPrimerFields const& fields, analysis::LegacyCriteria const& criteria,
                                std::chrono::year_month_day const& date) -> std::string {
    return legacy_multi_primer_report(analysis::with_extensions(fields), criteria, date);
}

auto legacy_multi_primer_report(analysis::LegacyMultiPrimerFields const& fields,
                                analysis::LegacyCriteria const& criteria, std::chrono::year_month_day const& date)
    -> std::string {
    auto const analysis = analysis::legacy_multi_primer(fields, criteria);
    auto report = title_line("Primer Design Report:  Multiple Primers", date);
    report += blank_lines(2);
    for (auto index = std::size_t{0}; index < analysis::legacy_multi_primer_fields; ++index) {
        report += multi_primer_section(index, fields);
    }
    report += line(rule);
    report += blank_lines(2);
    report += line("    Primer Balance and Cautions:              -A-     -B-     -C-     -D-");
    report += legacy_line_end;
    using Column = analysis::LegacyMultiPrimerColumn;
    report += multi_balance_row("Length", analysis, [](Column const& column) { return std::to_string(column.length); });
    report += multi_balance_row("% GC", analysis,
                                [](Column const& column) { return std::to_string(column.actual.gc_percent); });
    report += multi_balance_row("Tm C", analysis,
                                [](Column const& column) { return std::to_string(column.actual.tm_celsius); });
    report += multi_balance_row("Runs of G/C at 3' end", analysis,
                                [](Column const& column) { return std::to_string(column.actual.gc_run_3prime); });
    report += multi_balance_row("Runs of bases (any)", analysis,
                                [](Column const& column) { return std::to_string(column.actual.longest_run); });
    report += multi_balance_row("Hairpins", analysis,
                                [](Column const& column) { return std::string{column.cautions.hairpin ? "!" : "OK"}; });
    report += blank_lines(2);
    report += line("   Primer Dimers: (number of adjacent bases in largest duplex)");
    report += legacy_line_end;
    report += line("              A     B     C     D                  A     B     C     D");
    report += legacy_line_end;
    for (auto row = std::size_t{0}; row < analysis::legacy_multi_primer_fields; ++row) {
        report += matrix_row(row, analysis);
    }
    return report;
}


auto legacy_find_pair_report(FindPairReport const& find_pair, analysis::LegacyCriteria const& criteria,
                             std::chrono::year_month_day const& date) -> std::string {
    auto const& settings = find_pair.settings;
    auto report = title_line("Primer Design Report:  Find Pair", date);
    report += blank_lines(2);
    report += line(std::format("Source molecule:  {:<{}}{:>{}} bps", find_pair.molecule_name, molecule_name_width,
                               find_pair.molecule_size, molecule_size_width));
    report += blank_lines(2);
    report += line("Search parameters:");
    report += legacy_line_end;
    // the offsets as searched: cut at the molecule's ends (target - 1
    // bases before target A, the rest of the molecule after target B)
    auto const start_offset = std::min(settings.start_offset_a, settings.target_a > 0 ? settings.target_a - 1 : 0);
    auto const end_offset = std::min(settings.end_offset_b,
                                     find_pair.molecule_size > settings.target_b
                                         ? find_pair.molecule_size - settings.target_b
                                         : 0);
    report += line(std::format("   Primer A (normal strand): Target bp:{:6};   Start: Target - {} bases",
                               settings.target_a, start_offset));
    report += line(std::format("            Length:{:4}      5' Extension: - none -", settings.length));
    report += line(std::format("   Primer B (compl. strand): Target bp:{:6};   End:   Target + {} bases",
                               settings.target_b, end_offset));
    report += line(std::format("            Length:{:4}      5' Extension: - none -", settings.length));
    report += legacy_line_end;
    report += line(std::format("   Primer Balance:           GC Match: +/- {} GC %     TM Match: +/- {} C",
                               settings.gc_tolerance_percent, settings.tm_tolerance_celsius));
    report += legacy_line_end;
    report += short_criteria(criteria);
    report += legacy_line_end;
    report += line("Primer pairs requested:");
    report += legacy_line_end;
    report += line(std::format("       Position          GC(%)      Tm({}C)     Max Dimers   FalsePriming",
                               cp437_degree_sign));
    report += line("     -A-       -B-       A / B      A / B        3'  Any      TmA/TmB");
    report += legacy_line_end;
    for (auto const& pair : find_pair.best_pairs) {
        report += find_pair_row(pair);
    }
    return report;
}


auto legacy_homology_report(HomologyReport const& homology, std::chrono::year_month_day const& date) -> std::string {
    assert(homology.molecule_size > 0);
    auto const uniqueness = homology.kind == analysis::HomologyKind::uniqueness;
    auto report = title_line(uniqueness ? "Uniqueness of Primer Report" : "False Priming Report", date);
    report += blank_lines(2);
    report += line(rule);
    report += legacy_line_end;
    report += line(std::format("Primer Sequence:  5'-{}-3'", letters(homology.primer)));
    report += line(std::format("                  Bp pos:{:6}", homology.position));
    report += legacy_line_end;
    report += line(std::format("Source molecule:  {:<{}}{:>{}} bps", homology.molecule_name, molecule_name_width,
                               homology.molecule_size, molecule_size_width));
    report += blank_lines(2);
    report += line(rule);
    report += blank_lines(2);
    report += line(uniqueness ? "Search criteria defined:" : "False priming site criteria defined:");
    report += legacy_line_end;
    report += line(uniqueness ? "     Regions of homology (>= 7 bases) between primer and template"
                              : "     3' end homologies (>= 6 bases) between primer and template");
    // an empty line when sites with Tm <= 0 are kept
    report += line(homology.tm_filter == analysis::HomologyTmFilter::exclude_at_most_zero
                       ? "     Sites with melting temperatures <= 0 C are not reported"
                       : "");
    report += blank_lines(3);
    report += line(uniqueness ? "Molecule map with primer location (P) and homologous regions marked"
                              : "Molecule map with primer location (P) and false priming sites marked");
    report += legacy_line_end;
    report += homology_map(homology);
    report += blank_lines(3);
    report += line(uniqueness ? "Regions of homology found:" : "False priming sites found:");
    report += legacy_line_end;
    report += line("             Homologous bases                   Bp Pos          Tm C");
    report += legacy_line_end;
    for (auto const& site : homology.sites) {
        report += homology_row(homology, site);
    }
    return report;
}


auto legacy_match_report(MatchReport const& match, std::chrono::year_month_day const& date) -> std::string {
    auto const& settings = match.settings;
    auto report = title_line("Primer Design Report:  Match Active Primer", date);
    report += blank_lines(2);
    report += line(std::format("Source molecule:  {:<{}}{:>{}} bps", match.molecule_name, molecule_name_width,
                               match.molecule_size, molecule_size_width));
    report += blank_lines(2);
    report += line("Search parameters:");
    // primers B on the top strand: no strand named (recorded 2026-10-20,
    // match_variations/TRUN_50_600_top_strand.PRN)
    auto const strand = settings.strand == sequence::Strand::complement ? std::string_view{" (Complementary strand)"}
                                                                         : std::string_view{};
    report += line(std::format("   Target bp pos:{:6}{}", settings.target, strand));
    auto const [start_offset, end_offset] = match_offsets(match);
    report += line(std::format("   Start: Target - {} bases;   End: Target + {} bases", start_offset, end_offset));
    report += line(std::format("   Length:{:4}                 5' Extension: - none -", settings.length));
    report += line(std::format("   GC Match:  +/- {:2} GC %      Tm Match:  +/- {:2} C",
                               search::detail::form_field(settings.gc_tolerance_percent),
                               search::detail::form_field(settings.tm_tolerance_celsius)));
    report += blank_lines(2);
    report += line("Primer A (Original) Characteristics:");
    report += line("                                       -- Self Dimers --        -- Runs --");
    report += line(std::format("  Position              GC%   Tm{}C     3'  any                  any  3'G/C",
                               cp437_degree_sign));
    // PRIMER.exe marks nothing on the row of a complement-strand primer A
    // (tests/legacy/complement/match_c400.PRN: GC 45 without '!')
    report += match_primer_a_row(match.primer_a_position, analysis::legacy_actual_values(match.primer_a),
                                 match.criteria,
                                 match.primer_a_strand == sequence::Strand::complement ? PrimerAFlags::hidden
                                                                                       : PrimerAFlags::shown);
    report += blank_lines(2);
    report += line("Primer B (Matching) List:");
    report += line("                                         ---  Dimers  ---");
    report += line("             bps from                  -Self-       -Inter-     -- Runs --");
    report += line(std::format("  Position   Target     GC%   Tm{}C     3'  any      3'  any     any  3'G/C",
                               cp437_degree_sign));
    report += legacy_line_end;
    for (auto const& candidate : match.list) {
        report += match_row(candidate);
    }
    return report;
}

}  // namespace primerpicker::report
