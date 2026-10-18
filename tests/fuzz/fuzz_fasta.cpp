// libFuzzer harness: FASTA reader (core/io/molecule_file.hpp)

#include "core/io/molecule_file.hpp"
#include "fuzz_input.hpp"

#include <cstddef>
#include <cstdint>

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto const input = primerpicker::fuzz::text(primerpicker::fuzz::bytes(data, size));
    auto const result = primerpicker::io::parse_fasta(input);
    if (result.molecule.has_value()) {
        primerpicker::fuzz::check(result.error == primerpicker::io::MoleculeError::none);
        primerpicker::fuzz::check(not result.molecule->sequence.empty());
        primerpicker::fuzz::check(primerpicker::fuzz::is_upper_case_letters(result.molecule->sequence));
    }
    return 0;
}
