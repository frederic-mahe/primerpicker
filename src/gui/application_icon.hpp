#pragma once

#include <QIcon>

namespace primerpicker::gui {

// The icon of the windows (packaging/icons/, built into the program):
// the arrows alone from 16 to 48 px, the words and arrows from 64 px
[[nodiscard]] auto application_icon() -> QIcon;

}  // namespace primerpicker::gui
