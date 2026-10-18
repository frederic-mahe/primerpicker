#pragma once

#include <cstdint>

namespace primerpicker::sequence {

// A circular molecule (a plasmid) has no ends: its last base is
// followed by its first
enum struct Topology : std::uint8_t { linear, circular };

}  // namespace primerpicker::sequence
