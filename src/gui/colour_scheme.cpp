#include "gui/colour_scheme.hpp"

#include <QGuiApplication>
#include <QString>
#include <QStyleHints>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

auto colour_scheme_name(ColourScheme const scheme) -> QString {
    switch (scheme) {
        case ColourScheme::system:
            return u"system"_s;
        case ColourScheme::light:
            return u"light"_s;
        case ColourScheme::dark:
            return u"dark"_s;
        default:
            return u"system"_s;  // unreachable: all enumerators are handled
    }
}

auto colour_scheme_from_name(QString const& name) -> ColourScheme {
    if (name == colour_scheme_name(ColourScheme::light)) {
        return ColourScheme::light;
    }
    if (name == colour_scheme_name(ColourScheme::dark)) {
        return ColourScheme::dark;
    }
    return ColourScheme::system;
}

void apply_colour_scheme(ColourScheme const scheme) {
    auto const qt_scheme = [scheme] {
        switch (scheme) {
            case ColourScheme::light:
                return Qt::ColorScheme::Light;
            case ColourScheme::dark:
                return Qt::ColorScheme::Dark;
            case ColourScheme::system:
            default:  // unreachable: all enumerators are handled
                return Qt::ColorScheme::Unknown;  // Unknown: the system's
        }
    }();
    QGuiApplication::styleHints()->setColorScheme(qt_scheme);
}

}  // namespace primerpicker::gui
