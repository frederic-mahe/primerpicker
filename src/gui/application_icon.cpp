#include "gui/application_icon.hpp"

#include <QIcon>
#include <QString>

#include <array>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

auto application_icon() -> QIcon {
    // the resources of src/gui/CMakeLists.txt; Qt picks the size
    static constexpr auto files = std::array{
        u":/icons/primerpicker-small-16.png", u":/icons/primerpicker-small-24.png",
        u":/icons/primerpicker-small-32.png", u":/icons/primerpicker-small-48.png",
        u":/icons/primerpicker-64.png",       u":/icons/primerpicker-128.png",
        u":/icons/primerpicker-256.png"};
    auto icon = QIcon{};
    for (auto const* const file : files) {
        icon.addFile(QString{file});
    }
    return icon;
}

}  // namespace primerpicker::gui
