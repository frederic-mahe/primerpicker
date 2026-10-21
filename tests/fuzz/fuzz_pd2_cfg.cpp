// libFuzzer harness: legacy settings files (pd2.cfg; core/io/pd2_cfg.hpp).
// A file that parses is written back (write_pd2_cfg): that file has the
// same size and parses, and writing it again gives the same bytes.

#include "core/io/pd2_cfg.hpp"
#include "fuzz_input.hpp"

#include <cstddef>
#include <cstdint>
#include <span>

namespace {

constexpr auto max_weight = 5;

constexpr auto is_weight(int const weight) noexcept -> bool {
    return weight >= 0 and weight <= max_weight;
}

}  // namespace

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto const input = std::as_bytes(primerpicker::fuzz::bytes(data, size));
    auto const result = primerpicker::io::parse_pd2_cfg(input);
    using primerpicker::fuzz::check;
    check(result.defaults.has_value() == (result.error == primerpicker::io::Pd2CfgError::none));
    if (not result.defaults.has_value()) {
        return 0;
    }
    // the weights are usable by the Find Pair score (0 to 5)
    auto const& weights = result.defaults->weights;
    check(is_weight(weights.gc_mismatch) and is_weight(weights.tm_mismatch) and is_weight(weights.dimer_3prime) and
          is_weight(weights.dimer_any) and is_weight(weights.false_priming));
    auto const written = primerpicker::io::write_pd2_cfg(*result.defaults);
    check(written.size() == input.size());
    auto const reread = primerpicker::io::parse_pd2_cfg(written);
    check(reread.defaults.has_value());
    if (reread.defaults.has_value()) {
        check(primerpicker::io::write_pd2_cfg(*reread.defaults) == written);
    }
    return 0;
}
