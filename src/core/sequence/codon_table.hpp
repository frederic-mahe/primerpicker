#pragma once

#include <array>
#include <string_view>

namespace primerpicker::sequence {

// One line of the codon table of Primer Designer 2.0 (Mutagenize, F2):
// an amino acid, its one-letter code, and its codons, the letters in
// parentheses being alternatives for the third base (for instance
// TT(AG),CT(ACGT): TTA, TTG, CTA, CTC, CTG, CTT)
struct CodonTableEntry {
    std::string_view name;
    char code;
    std::string_view codons;
};

// The table as PRIMER.exe lists it (tests/legacy/mutagenize/codon_table.tsv,
// recorded 2026-10-19): the 20 amino acids by name, stop codons left out;
// the standard genetic code (see translate in genetic_code.hpp)
constexpr auto legacy_codon_table = std::array{
    CodonTableEntry{.name = "Ala", .code = 'A', .codons = "GC(ACGT)"},
    CodonTableEntry{.name = "Arg", .code = 'R', .codons = "CG(ACGT),AG(AG)"},
    CodonTableEntry{.name = "Asn", .code = 'N', .codons = "AA(CT)"},
    CodonTableEntry{.name = "Asp", .code = 'D', .codons = "GA(CT)"},
    CodonTableEntry{.name = "Cys", .code = 'C', .codons = "TG(CT)"},
    CodonTableEntry{.name = "Gln", .code = 'Q', .codons = "CA(AG)"},
    CodonTableEntry{.name = "Glu", .code = 'E', .codons = "GA(AG)"},
    CodonTableEntry{.name = "Gly", .code = 'G', .codons = "GG(ACGT)"},
    CodonTableEntry{.name = "His", .code = 'H', .codons = "CA(CT)"},
    CodonTableEntry{.name = "Ile", .code = 'I', .codons = "AT(ACT)"},
    CodonTableEntry{.name = "Leu", .code = 'L', .codons = "TT(AG),CT(ACGT)"},
    CodonTableEntry{.name = "Lys", .code = 'K', .codons = "AA(AG)"},
    CodonTableEntry{.name = "Met", .code = 'M', .codons = "ATG"},
    CodonTableEntry{.name = "Phe", .code = 'F', .codons = "TT(CT)"},
    CodonTableEntry{.name = "Pro", .code = 'P', .codons = "CC(ACGT)"},
    CodonTableEntry{.name = "Ser", .code = 'S', .codons = "TC(ACGT),AG(CT)"},
    CodonTableEntry{.name = "Thr", .code = 'T', .codons = "AC(ACGT)"},
    CodonTableEntry{.name = "Trp", .code = 'W', .codons = "TGG"},
    CodonTableEntry{.name = "Tyr", .code = 'Y', .codons = "TA(CT)"},
    CodonTableEntry{.name = "Val", .code = 'V', .codons = "GT(ACGT)"},
};

}  // namespace primerpicker::sequence
