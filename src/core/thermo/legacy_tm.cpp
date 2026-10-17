#include "core/thermo/legacy_tm.hpp"

#include "core/thermo/units.hpp"

#include <cassert>
#include <cmath>
#include <optional>

namespace primerpicker::thermo {

auto legacy_displayed_tm(Celsius const temperature) noexcept -> std::optional<int> {
    assert(std::isfinite(temperature.degrees));
    auto const truncated = static_cast<int>(std::floor(temperature.degrees));
    if (truncated <= 0) {
        return std::nullopt;
    }
    return truncated;
}

}  // namespace primerpicker::thermo
