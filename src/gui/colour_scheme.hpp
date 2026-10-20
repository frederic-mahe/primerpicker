#pragma once

#include <QString>

#include <cstdint>

namespace primerpicker::gui {

// The colours of the windows: those of the operating system's setting,
// or always light, or always dark (proposal §3.2; needs Qt 6.8)
enum struct ColourScheme : std::uint8_t { system, light, dark };

// The name kept in the settings file
[[nodiscard]] auto colour_scheme_name(ColourScheme scheme) -> QString;

// The scheme of a name kept in the settings file; an unknown or missing
// name follows the system
[[nodiscard]] auto colour_scheme_from_name(QString const& name) -> ColourScheme;

// Asks Qt for the scheme; platforms without a dark style keep theirs
void apply_colour_scheme(ColourScheme scheme);

}  // namespace primerpicker::gui
