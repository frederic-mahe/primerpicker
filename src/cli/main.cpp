// primerpicker: command-line front-end. It parses arguments and prints
// results; every value is computed by the core library.

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

namespace {

using primerpicker::analysis::evaluate;
using primerpicker::analysis::legacy_actual_values;
using primerpicker::analysis::LegacyActualValues;
using primerpicker::analysis::LegacyCriteria;
using primerpicker::analysis::LegacyVerdict;
using primerpicker::sequence::Nucleotide;

constexpr auto program_name = std::string_view{"primerpicker"};
constexpr auto usage = std::string_view{
    "Usage: primerpicker analyze PRIMER\n"
    "       primerpicker --help | --version\n"
    "\n"
    "analyze  analyse one primer (5' to 3', A C G T) with the criteria and\n"
    "         methods of Primer Designer 2.0, as its Primer Analysis screen\n"};

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

auto analysis_report(std::span<Nucleotide const> const primer) -> std::string {
    auto const criteria = LegacyCriteria{};
    auto const actual = legacy_actual_values(primer);
    auto const verdict = evaluate(actual, criteria);
    auto report = std::format("Sequence:      5'-{}-3'\n", letters(primer));
    report += std::format("Length:        {}\n", primer.size());
    report += "5' extension:  --none--\n";
    report += std::format("Composition:   {}\n\n", composition(primer));
    report += criteria_table(actual, criteria, verdict);
    report += '\n';
    report += hairpin_line(actual);
    report += std::format("Meets criteria: {}\n", verdict_word(verdict.all()));
    return report;
}

auto analyze(std::string_view const text) -> int {
    if (auto const invalid = text.find_first_not_of("ACGTacgt"); invalid != std::string_view::npos) {
        return fail(std::format("'{}' at position {} is not a nucleotide (A, C, G or T)", text[invalid],
                                invalid + 1));
    }
    if (text.size() < min_primer_length or text.size() > max_primer_length) {
        return fail(std::format("a primer has {} to {} bases", min_primer_length, max_primer_length));
    }
    auto const primer = primerpicker::sequence::to_nucleotides(text);
    if (not primer.has_value()) {
        return fail("not a nucleotide sequence");
    }
    std::cout << analysis_report(*primer);
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
    if (command == "analyze" and arguments.size() == 3) {
        return analyze(arguments[2]);
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
