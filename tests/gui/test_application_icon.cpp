#include "gui/application_icon.hpp"

#include <QIcon>
#include <QImage>
#include <QList>
#include <QSize>
#include <QString>

#include <catch2/catch_test_macros.hpp>

#include <array>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::application_icon;

// the drawings of packaging/icons/ (made by make_icons.sh); compared
// premultiplied, as pixmaps keep them (the edges round otherwise)
auto icon_file(QString const& name) -> QImage {
    return QImage{QString::fromUtf8(PRIMERPICKER_ICONS_DIR) + u"/png/"_s + name}.convertToFormat(
        QImage::Format_ARGB32_Premultiplied);
}

auto icon_image(QIcon const& icon, int const pixels) -> QImage {
    return icon.pixmap(QSize{pixels, pixels}).toImage().convertToFormat(QImage::Format_ARGB32_Premultiplied);
}

}  // namespace

TEST_CASE("The application icon has every size of packaging/icons", "[gui][icon]") {
    auto const icon = application_icon();
    REQUIRE_FALSE(icon.isNull());
    auto const sizes = icon.availableSizes();
    for (auto const pixels : std::array{16, 24, 32, 48, 64, 128, 256}) {
        CAPTURE(pixels);
        CHECK(sizes.contains(QSize{pixels, pixels}));
    }
}

TEST_CASE("The application icon shows the arrows alone up to 48 px, and the words from 64 px", "[gui][icon]") {
    auto const icon = application_icon();
    CHECK(icon_image(icon, 16) == icon_file(u"primerpicker-small-16.png"_s));
    CHECK(icon_image(icon, 48) == icon_file(u"primerpicker-small-48.png"_s));
    CHECK(icon_image(icon, 64) == icon_file(u"primerpicker-64.png"_s));
    CHECK(icon_image(icon, 256) == icon_file(u"primerpicker-256.png"_s));
}
