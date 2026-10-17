#include "core/analysis/legacy_criteria.hpp"
#include "core/sequence/nucleotide.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <charconv>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace {

using primerpicker::analysis::gc_run_at_3prime_end;
using primerpicker::analysis::longest_run_before_3prime_end;
using primerpicker::analysis::longest_self_complementary_run;
using primerpicker::analysis::self_dimer_at_3prime_end;
using primerpicker::sequence::Nucleotide;
using primerpicker::sequence::to_nucleotides;

auto nucleotides(std::string_view const letters) -> std::vector<Nucleotide> {
    auto converted = to_nucleotides(letters);
    if (not converted.has_value()) {
        throw std::runtime_error{"not a nucleotide sequence: " + std::string{letters}};
    }
    return *converted;
}

}  // namespace

// Reference values: PRIMER.exe, Primer/Create -> Direct Entry ->
// Inspect, "Actual" column, built-in defaults
// (tests/legacy/inspect/inspect.tsv, recorded 2026-10-17)

TEST_CASE("runs ignore the run at the 3' end", "[analysis][legacy]") {
    CHECK(longest_run_before_3prime_end(nucleotides("AAAATATTCT")) == 4);    // PRIMER.exe
    CHECK(longest_run_before_3prime_end(nucleotides("CAACTGCTATTT")) == 2);  // PRIMER.exe
    CHECK(longest_run_before_3prime_end(nucleotides("CCCGTCCC")) == 3);      // PRIMER.exe
    CHECK(longest_run_before_3prime_end(nucleotides("TTTTTTTTTTTTTTTTTTTTTTTT")) == 0);  // PRIMER.exe
}

TEST_CASE("G or C run at the 3' end", "[analysis][legacy]") {
    CHECK(gc_run_at_3prime_end(nucleotides("GGCCTCCGCC")) == 5);    // PRIMER.exe
    CHECK(gc_run_at_3prime_end(nucleotides("CAACTGCTATTT")) == 0);  // PRIMER.exe
    CHECK(gc_run_at_3prime_end(nucleotides("GCGCGCGCGCGC")) == 12); // PRIMER.exe
}

TEST_CASE("longest self-complementary run, any alignment", "[analysis][legacy]") {
    CHECK(longest_self_complementary_run(nucleotides("GGCCTCCGCC")) == 4);       // PRIMER.exe
    CHECK(longest_self_complementary_run(nucleotides("AAAATATTCT")) == 6);       // PRIMER.exe
    CHECK(longest_self_complementary_run(nucleotides("GTATTACATATTAAC")) == 4);  // PRIMER.exe
    CHECK(longest_self_complementary_run(nucleotides("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")) == 30);  // PRIMER.exe
}

TEST_CASE("self-dimer at the 3' end needs a partner overhang", "[analysis][legacy]") {
    CHECK(self_dimer_at_3prime_end(nucleotides("GGCCTCCGCC")) == 1);       // PRIMER.exe
    CHECK(self_dimer_at_3prime_end(nucleotides("GTATTACATATTAAC")) == 0);  // PRIMER.exe
    CHECK(self_dimer_at_3prime_end(nucleotides("CAACTGCTATTT")) == 2);     // PRIMER.exe
    CHECK(self_dimer_at_3prime_end(nucleotides("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC")) == 28);  // PRIMER.exe
}

namespace {

// columns: primer, dimers_3prime, dimers_any, runs, gc_run_3prime, ...
struct Recording {
    std::string primer;
    std::size_t dimers_3prime;
    std::size_t dimers_any;
    std::size_t runs;
    std::size_t gc_run_3prime;

    friend auto operator==(Recording const&, Recording const&) -> bool = default;
};

auto computed(std::string const& primer) -> Recording {
    auto const converted = nucleotides(primer);
    return {.primer = primer,
            .dimers_3prime = self_dimer_at_3prime_end(converted),
            .dimers_any = longest_self_complementary_run(converted),
            .runs = longest_run_before_3prime_end(converted),
            .gc_run_3prime = gc_run_at_3prime_end(converted)};
}

auto to_size(std::string const& text) -> std::size_t {
    auto value = std::size_t{0};
    auto const* const last = std::next(text.data(), std::ssize(text));
    auto const [end, error] = std::from_chars(text.data(), last, value);
    if (error != std::errc{} or end != last) {
        throw std::runtime_error{"unexpected number: " + text};
    }
    return value;
}

auto read_recordings(std::filesystem::path const& path) -> std::vector<Recording> {
    constexpr auto first_columns = 5U;
    auto recordings = std::vector<Recording>{};
    auto file = std::ifstream{path};
    auto line = std::string{};
    while (std::getline(file, line)) {
        if (line.starts_with('#')) {
            continue;
        }
        auto fields = std::vector<std::string>{};
        auto stream = std::istringstream{line};
        for (auto field = std::string{}; fields.size() < first_columns and std::getline(stream, field, '\t');) {
            fields.push_back(field);
        }
        if (fields.size() != first_columns) {
            throw std::runtime_error{"unexpected line: " + line};
        }
        recordings.push_back({.primer = fields[0],
                              .dimers_3prime = to_size(fields[1]),
                              .dimers_any = to_size(fields[2]),
                              .runs = to_size(fields[3]),
                              .gc_run_3prime = to_size(fields[4])});
    }
    return recordings;
}

}  // namespace

TEST_CASE("legacy criteria match every recorded primer", "[analysis][legacy]") {
    auto const recordings =
        read_recordings(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "inspect" / "inspect.tsv");
    REQUIRE(not recordings.empty());
    for (auto const& recording : recordings) {
        INFO("primer " << recording.primer);
        CHECK(computed(recording.primer) == recording);
    }
}

TEST_CASE("legacy criteria work at compile time", "[analysis][legacy][constexpr]") {
    // GGCCTCCGCC (PRIMER.exe, inspect.tsv): 3' dimers 1, any 4, runs 2, 3' GC run 5
    constexpr auto primer = std::array{Nucleotide::guanine, Nucleotide::guanine, Nucleotide::cytosine,
                                       Nucleotide::cytosine, Nucleotide::thymine, Nucleotide::cytosine,
                                       Nucleotide::cytosine, Nucleotide::guanine, Nucleotide::cytosine,
                                       Nucleotide::cytosine};
    STATIC_CHECK(self_dimer_at_3prime_end(primer) == 1);
    STATIC_CHECK(longest_self_complementary_run(primer) == 4);
    STATIC_CHECK(longest_run_before_3prime_end(primer) == 2);
    STATIC_CHECK(gc_run_at_3prime_end(primer) == 5);
}
