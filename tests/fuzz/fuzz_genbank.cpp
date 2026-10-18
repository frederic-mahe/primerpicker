// libFuzzer harness: GenBank reader (core/io/genbank.hpp)

#include "core/io/genbank.hpp"
#include "fuzz_input.hpp"

#include <cstddef>
#include <cstdint>

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto const input = primerpicker::fuzz::text(primerpicker::fuzz::bytes(data, size));
    auto const result = primerpicker::io::parse_genbank(input);
    if (result.molecule.has_value()) {
        primerpicker::fuzz::check(result.error == primerpicker::io::GenbankError::none);
        primerpicker::fuzz::check(primerpicker::fuzz::is_upper_case_letters(result.molecule->sequence));
        primerpicker::fuzz::check(result.molecule->sequence.size() <= input.size());
    }
    return 0;
}
