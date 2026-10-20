#include "core/report/legacy_report.hpp"

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_duplex.hpp"
#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_pair.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <format>
#include <iterator>
#include <span>
#include <string>
#include <string_view>

namespace primerpicker::report {

namespace {

using analysis::PrimerPair;
using sequence::Nucleotide;

constexpr auto title_width = std::size_t{54};  // "Report Date:" starts in this column
constexpr auto lines_after_title = std::size_t{5};
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

auto letters(std::span<Nucleotide const> const primer) -> std::string {
    auto text = std::string{};
    std::ranges::transform(primer, std::back_inserter(text), sequence::to_letter);
    return text;
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
    assert(left.size() + criterion.actual.size() <= criteria_actual_end);
    return line(std::format("{}{:>{}}{}{}", left, criterion.actual, criteria_actual_end - left.size(), verdict_gap,
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

auto duplex_lines(PrimerPair const& pair_of_primers, analysis::DuplexPlacement const& placement) -> std::string {
    auto const marks_column = duplex_column + placement.marks_offset;
    assert(marks_column >= 0);
    auto text = line(std::string{duplex_indent} + std::string{five_prime_label} + letters(pair_of_primers.first) +
                     std::string{three_prime_end});
    text += line(std::string(static_cast<std::size_t>(marks_column), ' ') + std::string(placement.count, ':'));
    text += line(second_primer_line(pair_of_primers.second, duplex_column + placement.second_3prime_offset));
    return text;
}


// Primer Pair report: values right-aligned on these columns (inclusive),
// each followed by its flag ('!' or a space)
constexpr auto pair_label_column = std::size_t{14};
constexpr auto pair_first_value_end = std::size_t{40};  // one past the value of A
constexpr auto pair_value_step = std::size_t{9};       // to the next primer
constexpr auto pair_comment_gap = std::string_view{"     "};
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
        assert(row.size() + value.text.size() <= end);
        row += std::format("{:>{}}{}", value.text, end - row.size(), flag(value.flagged));
        end += pair_value_step;
    }
    return row;
}

auto primer_header(char const name) -> std::string {
    if (name == 'A') {
        return "-----  P R I M E R   A  ---------------------------------------------------";
    }
    return "-----------------------------------------------  P R I M E R   B  ---------";
}

auto primer_summary(char const name, std::span<Nucleotide const> const primer, analysis::LegacyCriteria const& criteria)
    -> std::string {
    auto const actual = analysis::legacy_actual_values(primer);
    auto const verdict = analysis::evaluate(actual, criteria);
    auto const count = [primer](Nucleotide const nucleotide) { return std::ranges::count(primer, nucleotide); };
    auto text = line(primer_header(name));
    text += legacy_line_end;
    text += line(std::format("Sequence:     5'-{}-3'", letters(primer)));
    text += legacy_line_end;
    text += line(std::format("              Length: {:3};", primer.size()));
    text += line(std::format("              % GC: {:2};     Tm C: {:2};    Composition:  A:{}  C:{}  G:{}  T:{}",
                             actual.gc_percent, actual.tm_celsius, count(Nucleotide::adenine),
                             count(Nucleotide::cytosine), count(Nucleotide::guanine), count(Nucleotide::thymine)));
    text += legacy_line_end;
    text += line(std::format("Criteria:     % GC:  Min {:3}, Max {:3};    Tm C:  Min {:3}, Max {:3}",
                             criteria.gc_min_percent, criteria.gc_max_percent, criteria.tm_min_celsius,
                             criteria.tm_max_celsius));
    text += line(std::format("              Hairpin cutoff: {:.1f} kcal;", criteria.hairpin_cutoff.kcal_per_mol));
    text += line(std::format("              Reject >= {} matches at 3' end;  >= {} adj homol bases",
                             criteria.reject_dimer_3prime_matches, criteria.reject_dimer_adjacent_bases));
    text += line(std::format("              Reject >= {} base runs;  >= {} G/C at 3' end", criteria.reject_run_length,
                             criteria.reject_gc_run_3prime));
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

auto balance_row(Balance const& balance) -> std::string {
    auto const values = std::array{FlaggedValue{.text = std::to_string(balance.value_a), .flagged = false},
                                   FlaggedValue{.text = std::to_string(balance.value_b), .flagged = false}};
    return line(std::format("{}{}{}{} difference", pair_row(balance.label, values), pair_comment_gap,
                            balance.difference, balance.unit));
}

auto caution_row(std::string_view const label, FlaggedValue const& first, FlaggedValue const& second) -> std::string {
    return line(pair_row(label, std::array{first, second}));
}

// "OK" or "!" for each primer, without flags (nothing after the last)
auto hairpin_row(bool const hairpin_a, bool const hairpin_b) -> std::string {
    auto const text = [](bool const hairpin) { return hairpin ? "!" : "OK"; };
    auto const row = std::format("{:{}}{}", "", pair_label_column, "Hairpins");
    return line(std::format("{}{:>{}}{:>{}}", row, text(hairpin_a), pair_first_value_end - row.size(), text(hairpin_b),
                            pair_value_step));
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

auto multi_primer_section(std::size_t const index, std::span<Nucleotide const> const primer) -> std::string {
    auto text = line(std::format("-----  P R I M E R   {}  ---------------------------------------------------",
                                 multi_names.at(index)));
    text += legacy_line_end;
    text += line(std::format("Sequence:     5'-{}-3'", letters(primer)));
    text += legacy_line_end;
    text += line("              Source molecule:  none");
    text += line(std::format("              Length: {:3};", primer.size()));
    text += legacy_line_end;
    return text;
}

// A row of "Primer Balance and Cautions": one value per primer, no flag
template <typename Value>
auto multi_balance_row(std::string_view const label, analysis::LegacyMultiPrimer const& analysis, Value const value)
    -> std::string {
    auto canvas = Canvas{};
    canvas.put(multi_label_column, label);
    for (auto index = std::size_t{0}; index < analysis::legacy_multi_primer_fields; ++index) {
        auto const& column = analysis.columns.at(index);
        if (column.has_value()) {
            canvas.put_right(multi_first_value_last + (index * multi_value_step), value(*column));
        }
    }
    return canvas.line();
}

// Row ROW of both dimer matrices: dots below the diagonal, values with
// their '!' flag on and above it
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
        if (row > 0) {
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

}  // namespace

auto legacy_inspect_report(analysis::LegacyPrimer const& primer, analysis::LegacyCriteria const& criteria,
                           std::chrono::year_month_day const& date) -> std::string {
    auto const actual = analysis::legacy_actual_values(primer);
    auto const verdict = analysis::evaluate(actual, criteria);
    auto const& sequence = primer.sequence;
    auto const extension = sequence.first(primer.extension_length);
    auto const self = PrimerPair{.first = sequence, .second = sequence};
    auto const count = [sequence](Nucleotide const nucleotide) { return std::ranges::count(sequence, nucleotide); };

    auto report = title_line("Primer Design Report:  One Primer", date);
    report += blank_lines(lines_after_title);
    report += line("-----  P R I M E R   I N F O  --------------------------------------------");
    report += blank_lines(2);
    // the extension, then a dash, then the primer
    auto const shown = extension.empty() ? letters(sequence)
                                         : letters(extension) + "-" + letters(sequence.subspan(extension.size()));
    report += line(std::format("Sequence:     5'-{}-3'", shown));
    report += legacy_line_end;
    report += line(std::format("              Length: {:3};", sequence.size()));
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


auto legacy_pair_report(PrimerPair const& primers, analysis::LegacyCriteria const& criteria,
                        std::chrono::year_month_day const& date) -> std::string {
    auto const& [primer_a, primer_b] = primers;
    auto const actual_a = analysis::legacy_actual_values(primer_a);
    auto const actual_b = analysis::legacy_actual_values(primer_b);
    auto const cautions_a = analysis::legacy_pair_cautions(actual_a, criteria);
    auto const cautions_b = analysis::legacy_pair_cautions(actual_b, criteria);
    auto const values_ab = analysis::legacy_pair_values(primers);
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
    report += primer_summary('A', primer_a, criteria);
    report += primer_summary('B', primer_b, criteria);
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
    auto const analysis = analysis::legacy_multi_primer(fields, criteria);
    auto report = title_line("Primer Design Report:  Multiple Primers", date);
    report += blank_lines(2);
    for (auto index = std::size_t{0}; index < analysis::legacy_multi_primer_fields; ++index) {
        if (auto const& field = fields.at(index); field.has_value()) {
            report += multi_primer_section(index, *field);
        }
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

}  // namespace primerpicker::report
