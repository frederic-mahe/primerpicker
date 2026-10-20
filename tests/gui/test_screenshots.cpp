// Screenshots of the main window, for reviews and the README. Hidden
// test (not run by CTest): run it with
//   PRIMERPICKER_SCREENSHOT_DIR=DIR build/debug/tests/gui/gui_tests '[.screenshots]'

#include "gui/colour_scheme.hpp"
#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QApplication>
#include <QDir>
#include <QLineEdit>
#include <QString>
#include <QTest>
#include <QtEnvironmentVariables>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::ColourScheme;
using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

void save_analysis(QDir const& directory, ColourScheme const scheme) {
    primerpicker::gui::apply_colour_scheme(scheme);
    auto window = MainWindow{};
    constexpr auto width = 760;
    constexpr auto height = 560;
    window.resize(width, height);
    window.show();
    child<QAction>(window, u"direct_entry").trigger();
    auto& field = child<QLineEdit>(window, u"direct_entry_field");
    QTest::keyClicks(&field, u"CGTCCAGCACCTTGG"_s);
    QTest::keyClick(&field, Qt::Key_Return);
    QApplication::processEvents();
    auto const file = directory.filePath(u"analysis_%1.png"_s.arg(primerpicker::gui::colour_scheme_name(scheme)));
    CHECK(window.grab().save(file));
}

}  // namespace

TEST_CASE("screenshots of the Primer Analysis view", "[.screenshots]") {
    auto const directory = qEnvironmentVariable("PRIMERPICKER_SCREENSHOT_DIR");
    REQUIRE_FALSE(directory.isEmpty());
    // the offscreen platform has no dark palette; Fusion follows the
    // requested scheme on every platform
    QApplication::setStyle(u"Fusion"_s);
    save_analysis(QDir{directory}, ColourScheme::light);
    save_analysis(QDir{directory}, ColourScheme::dark);
    primerpicker::gui::apply_colour_scheme(ColourScheme::system);
}
