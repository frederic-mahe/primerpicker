#pragma once

#include <QObject>
#include <QString>
#include <QStringView>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_tostring.hpp>

#include <string>

// QString in Catch2's messages
template <>
struct Catch::StringMaker<QString> {
    static auto convert(QString const& text) -> std::string { return '"' + text.toStdString() + '"'; }
};

namespace primerpicker::gui::test {

// The descendant of PARENT named NAME (its objectName). A view, not a
// QString temporary: GCC warns that the returned reference might
// dangle.
template <typename Widget>
auto child(QObject const& parent, QStringView const name) -> Widget& {
    auto* const found = parent.findChild<Widget*>(name.toString());
    INFO(name.toString().toStdString());
    REQUIRE(found != nullptr);
    return *found;
}

}  // namespace primerpicker::gui::test
