// libFuzzer harness: any molecule file (GenBank, FASTA, plain sequence;
// core/io/molecule_file.hpp)

#include "core/io/molecule_file.hpp"
#include "core/sequence/nucleotide.hpp"
#include "fuzz_input.hpp"

#include <cstddef>
#include <cstdint>

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto const input = primerpicker::fuzz::text(primerpicker::fuzz::bytes(data, size));
    auto const result = primerpicker::io::parse_molecule_file({.text = input, .name = "FUZZ"});
    if (not result.molecule.has_value()) {
        primerpicker::fuzz::check(result.error != primerpicker::io::MoleculeError::none);
        return 0;
    }
    auto const& sequence = result.molecule->sequence;
    primerpicker::fuzz::check(not sequence.empty());
    primerpicker::fuzz::check(primerpicker::fuzz::is_upper_case_letters(sequence));
    // what the front-ends do next with a molecule
    auto const nucleotides = primerpicker::sequence::to_nucleotides(sequence);
    if (nucleotides.has_value()) {
        primerpicker::fuzz::check(nucleotides->size() == sequence.size());
    }
    return 0;
}
