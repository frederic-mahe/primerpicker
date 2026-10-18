// primerpicker: command-line front-end. It parses arguments and prints
// results; every value is computed by the core library.

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_pair.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <format>
#include <iostream>
#include <iterator>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace {

using primerpicker::analysis::evaluate;
using primerpicker::analysis::legacy_actual_values;
using primerpicker::analysis::LegacyActualValues;
using primerpicker::analysis::LegacyCriteria;
using primerpicker::analysis::LegacyPrimer;
using primerpicker::analysis::LegacyVerdict;
using primerpicker::analysis::PrimerPair;
using primerpicker::sequence::Nucleotide;

constexpr auto program_name = std::string_view{"primerpicker"};
constexpr auto usage = std::string_view{
    "Usage: primerpicker analyze [--extension EXTENSION] PRIMER\n"
    "       primerpicker pair PRIMER_A PRIMER_B\n"
    "       primerpicker --help | --version\n"
    "\n"
    "Primers are written 5' to 3' with A, C, G and T, and analysed with the\n"
    "criteria and methods of Primer Designer 2.0.\n"
    "analyze  one primer, with an optional 5' extension, as the Primer\n"
    "         Analysis screen\n"
    "pair     two primers, as the Primer Pair Analysis screen\n"};

// Primer Designer's Direct Entry field holds 50 bases
constexpr auto max_primer_length = std::size_t{50};
constexpr auto min_primer_length = std::size_t{2};

auto fail(std::string_view const message) -> int {
    std::cerr << program_name << ": " << message << '\n';
    return EXIT_FAILURE;
}

auto verdict_word(bool const meets) -> std::string_view {
    return meets ? "YES" : "NO";
}

auto letter(Nucleotide const nucleotide) -> char {
    switch (nucleotide) {
        case Nucleotide::adenine:
            return 'A';
        case Nucleotide::cytosine:
            return 'C';
        case Nucleotide::guanine:
            return 'G';
        case Nucleotide::thymine:
            return 'T';
        default:
            return '?';  // unreachable: all enumerators are handled
    }
}

auto letters(std::span<Nucleotide const> const primer) -> std::string {
    auto text = std::string{};
    std::ranges::transform(primer, std::back_inserter(text), letter);
    return text;
}

auto composition(std::span<Nucleotide const> const primer) -> std::string {
    auto const count = [primer](Nucleotide const nucleotide) { return std::ranges::count(primer, nucleotide); };
    return std::format("A:{}  C:{}  G:{}  T:{}", count(Nucleotide::adenine), count(Nucleotide::cytosine),
                       count(Nucleotide::guanine), count(Nucleotide::thymine));
}

auto criterion_line(std::string_view const name, std::string_view const setting, std::string_view const actual,
                    bool const meets) -> std::string {
    return std::format("{:<17}{:<41}{:>6}  {}\n", name, setting, actual, verdict_word(meets));
}

auto criteria_table(LegacyActualValues const& actual, LegacyCriteria const& criteria, LegacyVerdict const& verdict)
    -> std::string {
    auto table = std::format("{:<17}{:<41}{:>6}  {}\n", "Criteria", "Setting", "Actual", "Meets");
    table += criterion_line("% GC", std::format("min {}, max {}", criteria.gc_min_percent, criteria.gc_max_percent),
                            std::to_string(actual.gc_percent), verdict.gc);
    table += criterion_line("Tm °C", std::format("min {}, max {}", criteria.tm_min_celsius, criteria.tm_max_celsius),
                            std::to_string(actual.tm_celsius), verdict.tm);
    // as on the Primer Analysis screen, the hairpin energy is not in
    // this column (see the "Hairpin:" line)
    table += criterion_line("No hairpins", std::format("energy cutoff {:.1f} kcal", criteria.hairpin_cutoff.kcal_per_mol),
                            "--", verdict.no_hairpin);
    table += criterion_line("No 3' dimers",
                            std::format("reject >= {} matches at 3' end", criteria.reject_dimer_3prime_matches),
                            std::to_string(actual.dimer_3prime_matches), verdict.no_3prime_dimer);
    table += criterion_line("No dimers",
                            std::format("reject >= {} adjacent homologous bases", criteria.reject_dimer_adjacent_bases),
                            std::to_string(actual.dimer_adjacent_bases), verdict.no_dimer);
    table += criterion_line("No runs", std::format("reject >= {} base runs", criteria.reject_run_length),
                            std::to_string(actual.longest_run), verdict.no_runs);
    table += criterion_line("No 3' GC runs", std::format("reject >= {} G or C at 3' end", criteria.reject_gc_run_3prime),
                            std::to_string(actual.gc_run_3prime), verdict.no_3prime_gc_run);
    return table;
}

auto hairpin_line(LegacyActualValues const& actual) -> std::string {
    if (not actual.hairpin.has_value()) {
        return "Hairpin:       no hairpins found\n";
    }
    return std::format("Hairpin:       {:.1f} kcal\n", actual.hairpin->kcal_per_mol);
}

auto analysis_report(LegacyPrimer const& primer) -> std::string {
    auto const criteria = LegacyCriteria{};
    auto const actual = legacy_actual_values(primer);
    auto const verdict = evaluate(actual, criteria);
    auto const extension = primer.sequence.first(primer.extension_length);
    auto report = std::format("Sequence:      5'-{}-3'\n", letters(primer.sequence));
    report += std::format("Length:        {}\n", primer.sequence.size());
    report += std::format("5' extension:  {}\n", extension.empty() ? std::string{"--none--"} : letters(extension));
    report += std::format("Composition:   {}\n\n", composition(primer.sequence));
    report += criteria_table(actual, criteria, verdict);
    report += '\n';
    report += hairpin_line(actual);
    report += std::format("Meets criteria: {}\n", verdict_word(verdict.all()));
    return report;
}

// The nucleotides of a primer typed on the command line, or none after
// printing why on stderr
auto parse_sequence(std::string_view const text, std::string_view const what, std::size_t const min_length)
    -> std::optional<std::vector<Nucleotide>> {
    if (auto const invalid = text.find_first_not_of("ACGTacgt"); invalid != std::string_view::npos) {
        fail(std::format("{}: '{}' at position {} is not a nucleotide (A, C, G or T)", what, text[invalid],
                         invalid + 1));
        return std::nullopt;
    }
    if (text.size() < min_length or text.size() > max_primer_length) {
        fail(std::format("{}: a primer has {} to {} bases", what, min_length, max_primer_length));
        return std::nullopt;
    }
    return primerpicker::sequence::to_nucleotides(text);
}

struct AnalyzeArguments {
    std::string_view extension;
    std::string_view primer;
};

auto analyze(AnalyzeArguments const& arguments) -> int {
    auto const extension = parse_sequence(arguments.extension, "extension", 0);
    auto const primer = parse_sequence(arguments.primer, "primer", min_primer_length);
    if (not extension.has_value() or not primer.has_value()) {
        return EXIT_FAILURE;
    }
    auto sequence = *extension;
    sequence.insert(sequence.end(), primer->begin(), primer->end());
    std::cout << analysis_report(LegacyPrimer{.sequence = sequence, .extension_length = extension->size()});
    return EXIT_SUCCESS;
}

auto flagged(auto const value, bool const flag) -> std::string {
    return std::format("{}{}", value, flag ? "!" : "");
}

auto pair_line(std::string_view const name, std::string_view const first, std::string_view const second,
               std::string_view const third) -> std::string {
    auto line = std::format("{:<26}{:>6}  {:>6}  {}", name, first, second, third);
    while (line.ends_with(' ')) {
        line.pop_back();
    }
    return line + '\n';
}

auto difference_comment(int const difference, std::string_view const unit) -> std::string {
    return difference == 0 ? std::string{} : std::format("{}{} difference", difference, unit);
}

auto pair_report(std::span<Nucleotide const> const primer_a, std::span<Nucleotide const> const primer_b) -> std::string {
    auto const criteria = LegacyCriteria{};
    auto const actual_a = legacy_actual_values(primer_a);
    auto const actual_b = legacy_actual_values(primer_b);
    auto const cautions_a = primerpicker::analysis::legacy_pair_cautions(actual_a, criteria);
    auto const cautions_b = primerpicker::analysis::legacy_pair_cautions(actual_b, criteria);
    auto const pair = primerpicker::analysis::legacy_pair_values(PrimerPair{.first = primer_a, .second = primer_b});
    auto const dimer_flags = [&criteria](std::size_t const three_prime, std::size_t const adjacent) {
        return primerpicker::analysis::legacy_dimer_flags(three_prime, adjacent, criteria);
    };
    auto const flags_aa = dimer_flags(actual_a.dimer_3prime_matches, actual_a.dimer_adjacent_bases);
    auto const flags_bb = dimer_flags(actual_b.dimer_3prime_matches, actual_b.dimer_adjacent_bases);
    auto const flags_ab = dimer_flags(pair.dimer_3prime_matches, pair.dimer_adjacent_bases);
    auto const either = [](bool const lhs, bool const rhs, std::string_view const comment) {
        return lhs or rhs ? std::string{comment} : std::string{};
    };
    auto report = std::format("Primer A:      5'-{}-3'\n", letters(primer_a));
    report += std::format("Primer B:      5'-{}-3'\n\n", letters(primer_b));
    report += pair_line("Primer balance", "-A-", "-B-", "Comment");
    report += pair_line("Length", std::to_string(primer_a.size()), std::to_string(primer_b.size()),
                        difference_comment(pair.length_difference, " base"));
    report += pair_line("% GC", flagged(actual_a.gc_percent, cautions_a.gc), flagged(actual_b.gc_percent, cautions_b.gc),
                        difference_comment(pair.gc_difference_percent, " %"));
    report += pair_line("Tm °C", flagged(actual_a.tm_celsius, cautions_a.tm), flagged(actual_b.tm_celsius, cautions_b.tm),
                        difference_comment(pair.tm_difference_celsius, "°C"));
    report += "\nPrimer cautions\n";
    report += pair_line("Runs of G/C at 3' end", flagged(actual_a.gc_run_3prime, cautions_a.gc_run_3prime),
                        flagged(actual_b.gc_run_3prime, cautions_b.gc_run_3prime),
                        either(cautions_a.gc_run_3prime, cautions_b.gc_run_3prime, "Runs found"));
    report += pair_line("Runs of bases (any)", flagged(actual_a.longest_run, cautions_a.runs),
                        flagged(actual_b.longest_run, cautions_b.runs), either(cautions_a.runs, cautions_b.runs, "Runs found"));
    report += pair_line("Hairpins", cautions_a.hairpin ? "!" : "OK", cautions_b.hairpin ? "!" : "OK",
                        either(cautions_a.hairpin, cautions_b.hairpin, "Hairpin(s) found"));
    report += '\n';
    report += pair_line("Primer dimers", "A:A", "B:B", std::format("{:>6}", "A:B"));
    report += pair_line("3' end matches", flagged(actual_a.dimer_3prime_matches, flags_aa.three_prime),
                        flagged(actual_b.dimer_3prime_matches, flags_bb.three_prime),
                        std::format("{:>6}", flagged(pair.dimer_3prime_matches, flags_ab.three_prime)));
    report += pair_line("Any adj homol bases", flagged(actual_a.dimer_adjacent_bases, flags_aa.adjacent),
                        flagged(actual_b.dimer_adjacent_bases, flags_bb.adjacent),
                        std::format("{:>6}", flagged(pair.dimer_adjacent_bases, flags_ab.adjacent)));
    return report;
}

auto pair(std::string_view const text_a, std::string_view const text_b) -> int {
    auto const primer_a = parse_sequence(text_a, "primer A", min_primer_length);
    auto const primer_b = parse_sequence(text_b, "primer B", min_primer_length);
    if (not primer_a.has_value() or not primer_b.has_value()) {
        return EXIT_FAILURE;
    }
    std::cout << pair_report(*primer_a, *primer_b);
    return EXIT_SUCCESS;
}

auto run(std::span<char const* const> const arguments) -> int {
    if (arguments.size() < 2) {
        std::cerr << usage;
        return EXIT_FAILURE;
    }
    auto const command = std::string_view{arguments[1]};
    if (command == "--help") {
        std::cout << usage;
        return EXIT_SUCCESS;
    }
    if (command == "--version") {
        std::cout << program_name << ' ' << PRIMERPICKER_VERSION << '\n';
        return EXIT_SUCCESS;
    }
    // program name, command, then the command's arguments
    constexpr auto analyze_count = std::size_t{3};
    constexpr auto analyze_with_extension_count = std::size_t{5};
    constexpr auto pair_count = std::size_t{4};
    constexpr auto extension_option = std::string_view{"--extension"};
    if (command == "analyze" and arguments.size() == analyze_count) {
        return analyze({.extension = "", .primer = arguments[2]});
    }
    if (command == "analyze" and arguments.size() == analyze_with_extension_count and
        arguments[2] == extension_option) {
        return analyze({.extension = arguments[3], .primer = arguments[4]});
    }
    if (command == "pair" and arguments.size() == pair_count) {
        return pair(arguments[2], arguments[3]);
    }
    std::cerr << usage;
    return EXIT_FAILURE;
}

}  // namespace

auto main(int const argc, char const* const* const argv) -> int {
    try {
        return run(std::span{argv, static_cast<std::size_t>(argc)});
    } catch (std::exception const& error) {  // std::bad_alloc, std::format_error
        return fail(error.what());
    }
}
