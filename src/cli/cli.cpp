// primerpicker: command-line front-end. It parses arguments and prints
// results; every value is computed by the core library.

#include "cli/cli.hpp"

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_pair.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/enzymes/enzyme.hpp"
#include "core/enzymes/rebase_enzymes.hpp"
#include "core/enzymes/site_search.hpp"
#include "core/io/molecule_file.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/sequence/iupac.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "core/sequence/topology.hpp"

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
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
    "       primerpicker find-pair [--weights WEIGHTS] MOLECULE TARGET_A TARGET_B\n"
    "       primerpicker find-enzymes [--almost] MOLECULE ENZYME\n"
    "       primerpicker --help | --version\n"
    "\n"
    "Primers are written 5' to 3' with A, C, G and T, and analysed with the\n"
    "criteria and methods of Primer Designer 2.0.\n"
    "analyze  one primer, with an optional 5' extension, as the Primer\n"
    "         Analysis screen\n"
    "pair     two primers, as the Primer Pair Analysis screen\n"
    "find-pair\n"
    "         pairs of primers around two targets of a molecule file (GenBank,\n"
    "         FASTA or plain sequence), ranked as the Find Primer Pair list;\n"
    "         primer A is searched from TARGET_A - 100 to TARGET_A, primer B\n"
    "         (complement strand) from TARGET_B to TARGET_B + 100. WEIGHTS: five\n"
    "         digits from 0 to 5, for GC mismatch, Tm mismatch, 3' dimers, any\n"
    "         dimers and false priming (default: 11111)\n"
    "find-enzymes\n"
    "         sites of a restriction enzyme (REBASE name; old names such as\n"
    "         Asp718 or Sau3A are accepted) on both strands of a molecule file,\n"
    "         as the True Sites screen; with --almost, sites with at most one\n"
    "         mismatch, as the Almost Sites screen\n"};

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
    auto const values_ab = primerpicker::analysis::legacy_pair_values(PrimerPair{.first = primer_a, .second = primer_b});
    auto const dimer_flags = [&criteria](std::size_t const three_prime, std::size_t const adjacent) {
        return primerpicker::analysis::legacy_dimer_flags(three_prime, adjacent, criteria);
    };
    auto const flags_aa = dimer_flags(actual_a.dimer_3prime_matches, actual_a.dimer_adjacent_bases);
    auto const flags_bb = dimer_flags(actual_b.dimer_3prime_matches, actual_b.dimer_adjacent_bases);
    auto const flags_ab = dimer_flags(values_ab.dimer_3prime_matches, values_ab.dimer_adjacent_bases);
    auto const either = [](bool const lhs, bool const rhs, std::string_view const comment) {
        return lhs or rhs ? std::string{comment} : std::string{};
    };
    auto report = std::format("Primer A:      5'-{}-3'\n", letters(primer_a));
    report += std::format("Primer B:      5'-{}-3'\n\n", letters(primer_b));
    report += pair_line("Primer balance", "-A-", "-B-", "Comment");
    report += pair_line("Length", std::to_string(primer_a.size()), std::to_string(primer_b.size()),
                        difference_comment(values_ab.length_difference, " base"));
    report += pair_line("% GC", flagged(actual_a.gc_percent, cautions_a.gc), flagged(actual_b.gc_percent, cautions_b.gc),
                        difference_comment(values_ab.gc_difference_percent, " %"));
    report += pair_line("Tm °C", flagged(actual_a.tm_celsius, cautions_a.tm), flagged(actual_b.tm_celsius, cautions_b.tm),
                        difference_comment(values_ab.tm_difference_celsius, "°C"));
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
                        std::format("{:>6}", flagged(values_ab.dimer_3prime_matches, flags_ab.three_prime)));
    report += pair_line("Any adj homol bases", flagged(actual_a.dimer_adjacent_bases, flags_aa.adjacent),
                        flagged(actual_b.dimer_adjacent_bases, flags_bb.adjacent),
                        std::format("{:>6}", flagged(values_ab.dimer_adjacent_bases, flags_ab.adjacent)));
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

auto read_text_file(std::filesystem::path const& path) -> std::optional<std::string> {
    auto file = std::ifstream{path, std::ios::binary};
    if (not file) {
        return std::nullopt;
    }
    auto text = std::ostringstream{};
    text << file.rdbuf();
    return std::move(text).str();
}

struct LoadedMolecule {
    std::vector<Nucleotide> sequence;
    primerpicker::sequence::Topology topology;
};

// The nucleotides and topology of a molecule file, or none after
// printing why on stderr
auto read_molecule(std::filesystem::path const& path) -> std::optional<LoadedMolecule> {
    auto const text = read_text_file(path);
    if (not text.has_value()) {
        fail(std::format("cannot read {}", path.string()));
        return std::nullopt;
    }
    auto const name = path.stem().string();
    auto const parsed = primerpicker::io::parse_molecule_file({.text = *text, .name = name});
    auto sequence = parsed.molecule.has_value() ? primerpicker::sequence::to_nucleotides(parsed.molecule->sequence)
                                                : std::nullopt;
    if (not sequence.has_value()) {
        fail(std::format("{}: not a DNA molecule (GenBank, FASTA or plain sequence of A, C, G and T)", path.string()));
        return std::nullopt;
    }
    return LoadedMolecule{.sequence = std::move(*sequence), .topology = parsed.molecule->topology};
}

// A position typed on the command line (1 or more), or none after
// printing why on stderr
auto parse_position(std::string_view const text, std::string_view const what) -> std::optional<std::size_t> {
    auto const digits = std::string{text};
    auto position = std::size_t{0};
    auto const* const end = std::next(digits.data(), std::ssize(digits));
    auto const [stop, error] = std::from_chars(digits.data(), end, position);
    if (error != std::errc{} or stop != end or position == 0) {
        fail(std::format("{}: '{}' is not a position (1 or more)", what, text));
        return std::nullopt;
    }
    return position;
}

// Five digits from 0 to 5, or none after printing why on stderr
auto parse_weights(std::string_view const text) -> std::optional<primerpicker::search::FindPairWeights> {
    constexpr auto weight_count = std::size_t{5};
    constexpr auto max_weight = '5';
    auto const is_weight = [](char const digit) { return digit >= '0' and digit <= max_weight; };
    if (text.size() != weight_count or not std::ranges::all_of(text, is_weight)) {
        fail(std::format("weights: '{}' is not five digits from 0 to 5", text));
        return std::nullopt;
    }
    auto const weight = [text](std::size_t const index) { return text[index] - '0'; };
    return primerpicker::search::FindPairWeights{.gc_mismatch = weight(0),
                                                 .tm_mismatch = weight(1),
                                                 .dimer_3prime = weight(2),
                                                 .dimer_any = weight(3),
                                                 .false_priming = weight(4)};
}

auto site_tm(std::optional<int> const tm_celsius) -> std::string {
    return tm_celsius.has_value() ? std::to_string(*tm_celsius) : std::string{"-"};
}

auto find_pair_line(primerpicker::search::PairCandidate const& pair) -> std::string {
    return std::format("{:>8}{:>10}{:>8} / {:<2}{:>6} / {:<2}{:>8}{:>5}{:>9} / {}\n", pair.position_a, pair.position_b,
                       pair.primer_a.gc_percent, pair.primer_b.gc_percent, pair.primer_a.tm_celsius,
                       pair.primer_b.tm_celsius, pair.max_dimer_3prime_matches, pair.max_dimer_adjacent_bases,
                       site_tm(pair.false_priming_tm_a), site_tm(pair.false_priming_tm_b));
}

// the messages of Primer Designer 2.0 (tests/legacy/README.md, "Limits
// on primers and pairs (recorded 2026-10-19)")
constexpr auto too_many_primers_message =
    std::string_view{"More than 100 sites found.  Suggest you redefine matching limits.\n"};
constexpr auto pairs_capped_message = std::string_view{"More than 1000 possible pairs found!\n"};

auto find_pair_report(primerpicker::search::FindPairResult const& result,
                      primerpicker::search::FindPairWeights const& weights) -> std::string {
    if (result.limit == primerpicker::search::FindPairLimit::too_many_primers) {
        return std::string{too_many_primers_message};
    }
    auto report = std::format("Primer A: {} possible primers\n", result.possible_a);
    report += std::format("Primer B: {} possible primers\n", result.possible_b);
    report += std::format("Pairs:    {} possible\n", result.pairs.size());
    if (result.limit == primerpicker::search::FindPairLimit::pairs_capped) {
        report += pairs_capped_message;
    }
    report += '\n';
    report += "   Position          GC(%)      Tm(°C)     Max Dimers   FalsePriming\n";
    report += "     -A-       -B-       A / B      A / B        3'  Any      TmA/TmB\n";
    auto const ranked = primerpicker::search::legacy_rank_pairs(result.pairs, weights);
    return std::accumulate(ranked.begin(), ranked.end(), std::move(report),
                           [](std::string text, primerpicker::search::PairCandidate const& candidate) {
                               return std::move(text) + find_pair_line(candidate);
                           });
}

struct FindPairArguments {
    std::string_view weights;
    std::string_view molecule;
    std::string_view target_a;
    std::string_view target_b;
};

auto find_pair(FindPairArguments const& arguments) -> int {
    auto const weights = parse_weights(arguments.weights);
    auto const target_a = parse_position(arguments.target_a, "target A");
    auto const target_b = parse_position(arguments.target_b, "target B");
    if (not weights.has_value() or not target_a.has_value() or not target_b.has_value()) {
        return EXIT_FAILURE;
    }
    auto const molecule = read_molecule(std::filesystem::path{arguments.molecule});
    if (not molecule.has_value()) {
        return EXIT_FAILURE;
    }
    auto const result = primerpicker::search::legacy_find_pair(
        molecule->sequence, {.target_a = *target_a, .target_b = *target_b}, LegacyCriteria{});
    std::cout << find_pair_report(result, *weights);
    return EXIT_SUCCESS;
}

// "....A": a dot per matching base, the molecule's base where it differs
// (across the origin of a circular molecule)
auto match_text(primerpicker::enzymes::Site const& site, std::span<Nucleotide const> const molecule,
                std::size_t const length) -> std::string {
    auto text = std::string(length, '.');
    if (site.mismatch_offset.has_value()) {
        auto const index = (site.position - 1 + *site.mismatch_offset) % molecule.size();
        text.at(*site.mismatch_offset) = letter(molecule[index]);
    }
    return text;
}

auto site_line(primerpicker::enzymes::Site const& site, std::span<Nucleotide const> const molecule,
               std::size_t const length, primerpicker::enzymes::Mismatches const mismatches) -> std::string {
    auto const complement = site.strand == primerpicker::sequence::Strand::complement;
    if (mismatches == primerpicker::enzymes::Mismatches::none) {
        return std::format("{:>8}{}\n", site.position, complement ? "C" : "");
    }
    auto const text = match_text(site, molecule, length);
    return complement ? std::format("{:>8} ({})\n", site.position, text)
                      : std::format("{:>8}  {}\n", site.position, text);
}

auto enzyme_report(primerpicker::enzymes::Enzyme const& enzyme, LoadedMolecule const& loaded,
                   primerpicker::enzymes::Mismatches const mismatches) -> std::optional<std::string> {
    auto const molecule = std::span<Nucleotide const>{loaded.sequence};
    auto const pattern = primerpicker::enzymes::parse_site(enzyme.site);
    if (not pattern.has_value()) {
        return std::nullopt;
    }
    auto const sites = primerpicker::enzymes::find_sites(molecule, *pattern, mismatches, loaded.topology);
    auto const reversed = primerpicker::sequence::reverse_complement_iupac(enzyme.site);
    // as Primer Designer, the almost sites header shows the reverse
    // complement of a non-palindromic site
    auto const shown_complement = mismatches == primerpicker::enzymes::Mismatches::at_most_one and
                                          reversed.has_value() and *reversed != enzyme.site
                                      ? std::format("  ({})", reversed.value_or(""))
                                      : std::string{};
    auto report = std::format("Search for: {:<8} {}{}    Number of Sites = {}\n", enzyme.name, enzyme.site,
                              shown_complement, sites.size());
    return std::accumulate(sites.begin(), sites.end(), std::move(report),
                           [&](std::string text, primerpicker::enzymes::Site const& site) {
                               return std::move(text) + site_line(site, molecule, enzyme.site.size(), mismatches);
                           });
}

struct FindEnzymesArguments {
    primerpicker::enzymes::Mismatches mismatches;
    std::string_view molecule;
    std::string_view enzyme;
};

auto enzyme_list(std::vector<primerpicker::enzymes::Enzyme> const& enzymes) -> std::string {
    auto list = std::string{};
    for (auto const& enzyme : enzymes) {
        list += (list.empty() ? "" : ", ") + std::string{enzyme.name};
    }
    return list.empty() ? std::string{"none"} : list;
}

// On stderr: what an old or unsold enzyme name stands for in REBASE,
// and the commercially available enzymes with the same site
auto note_enzyme(std::string_view const typed, primerpicker::enzymes::ResolvedEnzyme const& resolved) -> void {
    using primerpicker::enzymes::Availability;
    auto const& enzyme = resolved.enzyme;
    if (typed != enzyme.name) {
        std::cerr << std::format("{}: {} is {} in REBASE\n", program_name, typed, enzyme.name);
    }
    if (resolved.availability == Availability::commercial) {
        return;
    }
    std::cerr << std::format("{}: {} is not commercially available; isoschizomers: {}; neoschizomers (other cuts): "
                             "{}\n",
                             program_name, enzyme.name,
                             enzyme_list(primerpicker::enzymes::commercial_isoschizomers(enzyme)),
                             enzyme_list(primerpicker::enzymes::commercial_neoschizomers(enzyme)));
}

auto find_enzymes(FindEnzymesArguments const& arguments) -> int {
    auto const resolved = primerpicker::enzymes::resolve_enzyme(arguments.enzyme);
    if (not resolved.has_value()) {
        return fail(std::format("enzyme: '{}' is not a REBASE enzyme (names are case sensitive, for instance EcoRI)",
                                arguments.enzyme));
    }
    auto const molecule = read_molecule(std::filesystem::path{arguments.molecule});
    if (not molecule.has_value()) {
        return EXIT_FAILURE;
    }
    note_enzyme(arguments.enzyme, *resolved);
    auto const& enzyme = resolved->enzyme;
    auto const report = enzyme_report(enzyme, *molecule, arguments.mismatches);
    if (not report.has_value()) {
        return fail(std::format("enzyme {}: unexpected recognition sequence {}", enzyme.name, enzyme.site));
    }
    std::cout << *report;
    return EXIT_SUCCESS;
}

auto dispatch(std::span<char const* const> const arguments) -> int {
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
    constexpr auto find_pair_count = std::size_t{5};
    constexpr auto find_pair_with_weights_count = std::size_t{7};
    constexpr auto weights_option = std::string_view{"--weights"};
    constexpr auto default_weights = std::string_view{"11111"};
    if (command == "find-pair" and arguments.size() == find_pair_count) {
        return find_pair(
            {.weights = default_weights, .molecule = arguments[2], .target_a = arguments[3], .target_b = arguments[4]});
    }
    if (command == "find-pair" and arguments.size() == find_pair_with_weights_count and
        arguments[2] == weights_option) {
        auto const options = arguments.subspan(3);  // after "--weights"
        return find_pair({.weights = options[0], .molecule = options[1], .target_a = options[2], .target_b = options[3]});
    }
    constexpr auto find_enzymes_count = std::size_t{4};
    constexpr auto find_almost_enzymes_count = std::size_t{5};
    constexpr auto almost_option = std::string_view{"--almost"};
    if (command == "find-enzymes" and arguments.size() == find_enzymes_count) {
        return find_enzymes({.mismatches = primerpicker::enzymes::Mismatches::none,
                             .molecule = arguments[2],
                             .enzyme = arguments[3]});
    }
    if (command == "find-enzymes" and arguments.size() == find_almost_enzymes_count and
        arguments[2] == almost_option) {
        return find_enzymes({.mismatches = primerpicker::enzymes::Mismatches::at_most_one,
                             .molecule = arguments[3],
                             .enzyme = arguments[4]});
    }
    std::cerr << usage;
    return EXIT_FAILURE;
}

}  // namespace

namespace primerpicker::cli {

auto run(std::span<char const* const> const arguments) -> int {
    try {
        return dispatch(arguments);
    } catch (std::exception const& error) {  // std::bad_alloc, std::format_error
        return fail(error.what());
    }
}

}  // namespace primerpicker::cli
