#include "core/report/legacy_report.hpp"

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_duplex.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"

#include <algorithm>
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
    report += line(std::format("Sequence:     5'-{}-3'", letters(sequence)));
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

}  // namespace primerpicker::report
