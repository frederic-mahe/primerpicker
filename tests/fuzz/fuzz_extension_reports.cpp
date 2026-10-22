// libFuzzer harness: Primer Pair and Multi Primer analyses and reports
// of primers with a 5' extension (core/analysis/legacy_pair.hpp,
// legacy_multi_primer.hpp, core/report/legacy_report.hpp). The first
// byte says which of the fields A to D hold a primer (its low four
// bits); then, for each filled field, an extension byte (0 to
// max_extension_length bases), a length byte (2 to max_primer_length
// bases after the extension) and the bases (byte % 4). Missing bytes
// leave the remaining fields empty.

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_pair.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/report/legacy_report.hpp"
#include "core/sequence/nucleotide.hpp"
#include "fuzz_input.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace {

using primerpicker::analysis::legacy_multi_primer_fields;
using primerpicker::analysis::LegacyMultiPrimerFields;
using primerpicker::analysis::LegacyPrimer;
using primerpicker::sequence::Nucleotide;

constexpr auto min_primer_length = std::size_t{2};
constexpr auto max_primer_length = std::size_t{60};
constexpr auto max_extension_length = std::size_t{12};
constexpr auto bases = std::array{Nucleotide::adenine, Nucleotide::cytosine, Nucleotide::guanine, Nucleotide::thymine};
constexpr auto date = std::chrono::year_month_day{std::chrono::year{2026}, std::chrono::October, std::chrono::day{22}};

}  // namespace

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto input = primerpicker::fuzz::bytes(data, size);
    if (input.empty()) {
        return 0;
    }
    auto const filled = input.front();
    input = input.subspan(1);
    auto primers = std::array<std::vector<Nucleotide>, legacy_multi_primer_fields>{};
    auto fields = LegacyMultiPrimerFields{};
    constexpr auto field_header = std::size_t{2};  // extension and length bytes
    for (auto field = std::size_t{0}; field < legacy_multi_primer_fields; ++field) {
        if (((filled >> field) & 1U) == 0 or input.size() < field_header) {
            continue;
        }
        auto const extension = input[0] % (max_extension_length + 1);
        auto const length = extension + min_primer_length + (input[1] % (max_primer_length - min_primer_length + 1));
        input = input.subspan(field_header);
        if (input.size() < length) {
            break;
        }
        std::ranges::transform(input.first(length), std::back_inserter(primers.at(field)),
                               [](std::uint8_t const byte) { return bases.at(byte % bases.size()); });
        input = input.subspan(length);
        fields.at(field) = LegacyPrimer{.sequence = primers.at(field), .extension_length = extension};
    }

    using primerpicker::fuzz::check;
    auto const criteria = primerpicker::analysis::LegacyCriteria{};
    auto const analysis = primerpicker::analysis::legacy_multi_primer(fields, criteria);
    for (auto row = std::size_t{0}; row < legacy_multi_primer_fields; ++row) {
        auto const& field = fields.at(row);
        auto const& column = analysis.columns.at(row);
        check(column.has_value() == field.has_value());
        if (not field.has_value() or not column.has_value()) {
            continue;
        }
        // the length counts the extension; the Tm is that of the primer
        // without it
        check(column->length == field->sequence.size());
        auto const without = primerpicker::analysis::legacy_actual_values(field->sequence.subspan(field->extension_length));
        check(column->actual.tm_celsius == without.tm_celsius);
    }
    auto const multi_report = primerpicker::report::legacy_multi_primer_report(fields, criteria, date);
    check(primerpicker::fuzz::is_crlf_text(multi_report));

    // the Primer Pair of fields A and B
    auto const& field_a = fields.at(0);
    auto const& field_b = fields.at(1);
    if (field_a.has_value() and field_b.has_value()) {
        auto const pair = primerpicker::analysis::LegacyPrimerPair{.first = *field_a, .second = *field_b};
        auto const values = primerpicker::analysis::legacy_pair_values(pair);
        // the dimers of the whole primers
        auto const whole = primerpicker::analysis::legacy_pair_values(
            primerpicker::analysis::PrimerPair{.first = pair.first.sequence, .second = pair.second.sequence});
        check(values.dimer_3prime_matches == whole.dimer_3prime_matches);
        check(values.dimer_adjacent_bases == whole.dimer_adjacent_bases);
        auto const pair_report = primerpicker::report::legacy_pair_report(pair, criteria, date);
        check(primerpicker::fuzz::is_crlf_text(pair_report));
        // primer A's sequence is written with a dash after its extension
        constexpr auto sequence_mark = std::string_view{"Sequence:     5'-"};
        auto const first = pair_report.find(sequence_mark);
        check(first != std::string::npos);
        if (first != std::string::npos) {
            auto const start = first + sequence_mark.size();
            auto const end = pair_report.find("-3'", start);
            check(end != std::string::npos);
            auto const shown = pair_report.substr(start, end - start);
            check((shown.find('-') != std::string::npos) == (pair.first.extension_length > 0));
            check(shown.size() == pair.first.sequence.size() + (pair.first.extension_length > 0 ? 1 : 0));
        }
    }
    return 0;
}
