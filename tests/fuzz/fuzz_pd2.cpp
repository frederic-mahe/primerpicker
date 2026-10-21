// libFuzzer harness: legacy primer files (.PD2; core/io/pd2.hpp). A file
// that parses is written back (write_pd2): that file parses, and writing
// it again gives the same bytes.

#include "core/io/pd2.hpp"
#include "fuzz_input.hpp"

#include <cstddef>
#include <cstdint>
#include <span>

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto const input = std::as_bytes(primerpicker::fuzz::bytes(data, size));
    auto const result = primerpicker::io::parse_pd2(input);
    primerpicker::fuzz::check(result.primer.has_value() == (result.error == primerpicker::io::Pd2Error::none));
    if (not result.primer.has_value()) {
        return 0;
    }
    auto const written = primerpicker::io::write_pd2(*result.primer);
    auto const reread = primerpicker::io::parse_pd2(written);
    primerpicker::fuzz::check(reread.primer.has_value());
    if (reread.primer.has_value()) {
        primerpicker::fuzz::check(primerpicker::io::write_pd2(*reread.primer) == written);
    }
    return 0;
}
