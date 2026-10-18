// libFuzzer harness: recognition sequences and site search
// (core/enzymes/site_search.hpp, core/sequence/iupac.hpp). The input is
// a recognition sequence, a newline, then a molecule.

#include "core/enzymes/site_search.hpp"
#include "core/sequence/iupac.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/topology.hpp"
#include "fuzz_input.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string_view>

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    using primerpicker::enzymes::Mismatches;
    auto const input = primerpicker::fuzz::text(primerpicker::fuzz::bytes(data, size));
    auto const text = std::string_view{input};
    auto const newline = text.find('\n');
    auto const site = text.substr(0, newline);
    auto const molecule_text = newline == std::string_view::npos ? std::string_view{} : text.substr(newline + 1);

    auto const pattern = primerpicker::enzymes::parse_site(site);
    auto const reversed = primerpicker::sequence::reverse_complement_iupac(site);
    primerpicker::fuzz::check(reversed.has_value() == (pattern.has_value() or site.empty()));
    if (reversed.has_value()) {
        auto const twice = primerpicker::sequence::reverse_complement_iupac(*reversed);
        primerpicker::fuzz::check(twice.has_value() and twice->size() == site.size());
    }
    auto const molecule = primerpicker::sequence::to_nucleotides(molecule_text);
    if (not pattern.has_value() or not molecule.has_value()) {
        return 0;
    }
    using primerpicker::sequence::Topology;
    for (auto const topology : {Topology::linear, Topology::circular}) {
        auto const exact = primerpicker::enzymes::find_sites(*molecule, *pattern, Mismatches::none, topology);
        auto const almost = primerpicker::enzymes::find_sites(*molecule, *pattern, Mismatches::at_most_one, topology);
        primerpicker::fuzz::check(exact.size() <= almost.size());
        auto const last_start =
            topology == Topology::circular ? molecule->size() : molecule->size() + 1 - pattern->size();
        auto const fits = [&](primerpicker::enzymes::Site const& found) {
            return found.position >= 1 and found.position <= last_start and
                   (not found.mismatch_offset.has_value() or *found.mismatch_offset < pattern->size());
        };
        primerpicker::fuzz::check(std::ranges::all_of(almost, fits));
        primerpicker::fuzz::check(std::ranges::is_sorted(almost, {}, &primerpicker::enzymes::Site::position));
    }
    return 0;
}
