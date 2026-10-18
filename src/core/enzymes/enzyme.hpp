#pragma once

#include <array>
#include <string_view>

namespace primerpicker::enzymes {

// A cut, as REBASE writes it in its EMBOSS files (emboss_e.###): the
// cut is to the right of the given base on each strand, bases being
// numbered ... -2 -1 1 2 ... from the first base of the recognition
// sequence (1), on the top strand; 0 means no cut
struct Cut {
    int top_strand;
    int bottom_strand;

    friend constexpr auto operator==(Cut const&, Cut const&) -> bool = default;
};

// A restriction enzyme: its recognition sequence (IUPAC letters, 5' to
// 3') and where it cuts. cut_count is 2 (one cut on each strand), 4
// (two cuts on each strand) or 0 (unknown); unused cuts are {0, 0}.
// suppliers: REBASE supplier codes (data/rebase/emboss_s.###).
struct Enzyme {
    std::string_view name;
    std::string_view site;
    int cut_count;
    std::array<Cut, 2> cuts;
    std::string_view suppliers;
};

}  // namespace primerpicker::enzymes
