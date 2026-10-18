// libFuzzer harness: legacy primer files (.PD2; core/io/pd2.hpp)

#include "core/io/pd2.hpp"
#include "fuzz_input.hpp"

#include <cstddef>
#include <cstdint>
#include <span>

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto const input = std::as_bytes(primerpicker::fuzz::bytes(data, size));
    auto const result = primerpicker::io::parse_pd2(input);
    primerpicker::fuzz::check(result.primer.has_value() == (result.error == primerpicker::io::Pd2Error::none));
    return 0;
}
