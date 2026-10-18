#pragma once

#include "core/sequence/nucleotide.hpp"
#include "core/thermo/legacy_tm.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace primerpicker::analysis {

enum struct Strand : std::uint8_t { top, complement };

// Uniqueness: stretches of 7 bases or more anywhere in the primer;
// false priming: stretches of 6 bases or more at the primer's 3' end
enum struct HomologyKind : std::uint8_t { uniqueness, false_priming };

struct HomologySite {
    std::size_t primer_start;  // 0-based start of the matching stretch in the primer
    std::size_t length;        // bases
    std::size_t position;      // 1-based molecule position facing the stretch's 5'-most base
    Strand strand;             // complement: the match is on the complement strand
    int tm_celsius;            // legacy Tm of the stretch, truncated
};

namespace detail {

constexpr auto min_uniqueness_length = std::size_t{7};
constexpr auto min_false_priming_length = std::size_t{6};

// The molecule read 5' to 3' on one strand
class StrandView {
public:
    constexpr StrandView(std::span<sequence::Nucleotide const> const molecule, Strand const strand) noexcept
        : molecule_{molecule}, strand_{strand} {}

    [[nodiscard]] constexpr auto size() const noexcept -> std::size_t { return molecule_.size(); }

    [[nodiscard]] constexpr auto operator[](std::size_t const index) const noexcept -> sequence::Nucleotide {
        if (strand_ == Strand::top) {
            return molecule_[index];
        }
        return sequence::complement(molecule_[molecule_.size() - 1 - index]);
    }

    // 1-based top-strand coordinate of a base of this strand
    [[nodiscard]] constexpr auto position(std::size_t const index) const noexcept -> std::size_t {
        return strand_ == Strand::top ? index + 1 : molecule_.size() - index;
    }

    [[nodiscard]] constexpr auto strand() const noexcept -> Strand { return strand_; }

private:
    std::span<sequence::Nucleotide const> molecule_;
    Strand strand_;
};

// Truncated legacy Tm; not positive values are reported as 0
constexpr auto truncated_tm(std::span<sequence::Nucleotide const> const stretch) noexcept -> int {
    auto const stretch_tm = thermo::legacy_primer_tm(stretch).degrees;
    return stretch_tm > 0.0 ? static_cast<int>(stretch_tm) : 0;  // toward zero: the floor when positive
}

// Maximal matches of the primer on one strand, in primer then molecule
// order
constexpr auto append_sites(std::span<sequence::Nucleotide const> const primer, StrandView const& strand,
                            HomologyKind const kind, std::vector<HomologySite>& sites) -> void {
    auto const min_length = kind == HomologyKind::uniqueness ? min_uniqueness_length : min_false_priming_length;
    for (auto start = std::size_t{0}; start < primer.size(); ++start) {
        for (auto offset = std::size_t{0}; offset < strand.size(); ++offset) {
            auto const extends_a_match = start > 0 and offset > 0 and primer[start - 1] == strand[offset - 1];
            if (primer[start] != strand[offset] or extends_a_match) {
                continue;
            }
            auto length = std::size_t{0};
            while (start + length < primer.size() and offset + length < strand.size() and
                   primer[start + length] == strand[offset + length]) {
                ++length;
            }
            auto const reaches_3prime_end = start + length == primer.size();
            // a full-length match is the primer's own site
            if (length < min_length or length == primer.size() or
                (kind == HomologyKind::false_priming and not reaches_3prime_end)) {
                continue;
            }
            auto const site_tm = truncated_tm(primer.subspan(start, length));
            // sites with a Tm <= 0 °C are excluded (PRIMER.exe default)
            if (site_tm <= 0) {
                continue;
            }
            sites.push_back(HomologySite{.primer_start = start,
                                         .length = length,
                                         .position = strand.position(offset),
                                         .strand = strand.strand(),
                                         .tm_celsius = site_tm});
        }
    }
}

}  // namespace detail

// The Uniqueness or False Priming sites of Primer Designer 2.0, with
// its default settings (tests/legacy/README.md, "Uniqueness and false
// priming"): exact maximal matches on both strands, the primer's own
// site excluded, sites with a Tm <= 0 °C excluded, sorted by decreasing
// Tm (top-strand sites first among equal Tm). The primer must have at
// least two nucleotides. Not noexcept: the vector allocates.
// C++26 refactoring: make constexpr (std::ranges::stable_sort)
[[nodiscard]] inline auto legacy_homology_sites(std::span<sequence::Nucleotide const> const primer,
                                                   std::span<sequence::Nucleotide const> const molecule,
                                                   HomologyKind const kind) -> std::vector<HomologySite> {
    auto sites = std::vector<HomologySite>{};
    for (auto const strand : {Strand::top, Strand::complement}) {
        detail::append_sites(primer, detail::StrandView{molecule, strand}, kind, sites);
    }
    std::ranges::stable_sort(sites, [](HomologySite const& lhs, HomologySite const& rhs) {
        return lhs.tm_celsius > rhs.tm_celsius;
    });
    return sites;
}

}  // namespace primerpicker::analysis
