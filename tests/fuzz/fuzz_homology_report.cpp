// libFuzzer harness: the Uniqueness and False Priming search and reports
// (core/analysis/legacy_homology.hpp, core/report/legacy_report.cpp) on
// fuzzed molecules. The first bytes give the primer's position and
// length and the settings, each further byte a base (byte % 4), up to
// the largest molecule PRIMER.exe loads.

#include "core/analysis/legacy_homology.hpp"
#include "core/report/legacy_report.hpp"
#include "core/sequence/nucleotide.hpp"
#include "fuzz_input.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace {

using primerpicker::analysis::HomologyKind;
using primerpicker::analysis::HomologyOrder;
using primerpicker::analysis::HomologyTmFilter;
using primerpicker::sequence::Nucleotide;

constexpr auto setting_bytes = std::size_t{4};
constexpr auto max_molecule_length = std::size_t{54999};  // PRIMER.exe refuses 55,000 bases (README)
constexpr auto min_primer_length = std::size_t{2};
constexpr auto max_primer_length = std::size_t{50};
constexpr auto byte_bits = 8U;
constexpr auto map_scale_divisor = std::size_t{73};  // one character per L / 73 + 1 bases (README)
constexpr auto lines_from_map_title = 6;   // the title, an empty line, P>, the molecule, two strands
constexpr auto rows_after_map_title = 13;  // the map, the list's title and header
constexpr auto max_map_line = std::size_t{76};  // a mark in column 3 + (L - 1) / s, at most 75
constexpr auto bases = std::array{Nucleotide::adenine, Nucleotide::cytosine, Nucleotide::guanine, Nucleotide::thymine};

auto lines(std::string_view text) -> std::vector<std::string_view> {
    auto result = std::vector<std::string_view>{};
    for (auto end = text.find("\r\n"); end != std::string_view::npos; end = text.find("\r\n")) {
        result.push_back(text.substr(0, end));
        text.remove_prefix(end + 2);
    }
    return result;
}

}  // namespace

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto const input = primerpicker::fuzz::bytes(data, size);
    if (input.size() <= setting_bytes) {
        return 0;
    }
    auto const settings_input = input.first(setting_bytes);
    auto const molecule_input =
        input.subspan(setting_bytes).first(std::min(input.size() - setting_bytes, max_molecule_length));
    auto const length = min_primer_length + (settings_input[2] % (max_primer_length - min_primer_length + 1));
    if (molecule_input.size() < length) {
        return 0;
    }
    auto molecule = std::vector<Nucleotide>{};
    molecule.reserve(molecule_input.size());
    std::ranges::transform(molecule_input, std::back_inserter(molecule),
                           [](std::uint8_t const byte) { return bases.at(byte % bases.size()); });
    auto const value = static_cast<std::size_t>(settings_input[0]) | (static_cast<std::size_t>(settings_input[1]) << byte_bits);
    auto const position = 1 + (value % (molecule.size() - length + 1));
    auto const flags = settings_input[3];
    auto const kind = (flags & 1U) == 0 ? HomologyKind::uniqueness : HomologyKind::false_priming;
    auto const settings = primerpicker::analysis::HomologySettings{
        .tm_filter = (flags & 2U) == 0 ? HomologyTmFilter::exclude_at_most_zero : HomologyTmFilter::keep_all,
        .order = (flags & 4U) == 0 ? HomologyOrder::by_tm : HomologyOrder::by_position};

    auto const primer = std::span{molecule}.subspan(position - 1, length);
    auto const sites = primerpicker::analysis::legacy_homology_sites(primer, molecule, kind, settings);
    auto const report = primerpicker::report::legacy_homology_report(
        {.molecule_name = "FUZZ",
         .molecule_size = molecule.size(),
         .primer = primer,
         .position = position,
         .kind = kind,
         .tm_filter = settings.tm_filter,
         .sites = sites},
        std::chrono::year_month_day{std::chrono::year{2026}, std::chrono::October, std::chrono::day{20}});

    using primerpicker::fuzz::check;
    check(primerpicker::fuzz::is_crlf_text(report));
    // the map: P>, the molecule, the two strands, never wider than the page
    auto const all = lines(report);
    auto const map = std::ranges::find_if(all, [](std::string_view const line) { return line.starts_with("Molecule map"); });
    check(std::distance(map, all.end()) >= lines_from_map_title);
    auto const scale = (molecule.size() / map_scale_divisor) + 1;
    auto const map_lines = std::span{std::next(map, 2), 4};
    check(map_lines[1] == std::string(molecule.size() / scale, '='));
    check(std::ranges::all_of(map_lines, [](std::string_view const line) { return line.size() <= max_map_line; }));
    // one row per site, after the list's header
    check(static_cast<std::size_t>(std::distance(std::next(map, rows_after_map_title), all.end())) == sites.size());
    return 0;
}
