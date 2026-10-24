#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QDialog>
#include <QLabel>
#include <QPushButton>
#include <QString>
#include <QTabWidget>
#include <QTextBrowser>
#include <QtVersion>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

void open_about(MainWindow& window) {
    child<QAction>(window, u"about").trigger();
    REQUIRE(child<QDialog>(window, u"about_dialog").isVisible());
}

}  // namespace

TEST_CASE("Help -> About PrimerPicker names the program, its licence and the Qt version", "[gui][about]") {
    auto window = MainWindow{};
    open_about(window);
    auto const text = child<QLabel>(window, u"about_text").text();
    CHECK(text.contains(u"PrimerPicker"_s));
    CHECK(text.contains(u"GNU General Public License"_s));
    CHECK(text.contains(QString::fromLatin1(qVersion())));
    CHECK(text.contains(u"LGPL"_s));
    CHECK(text.contains(u"REBASE"_s));
    CHECK(child<QPushButton>(window, u"about_qt_button").isVisibleTo(&window));
}

TEST_CASE("the About dialog holds the licences of the packages", "[gui][about]") {
    auto window = MainWindow{};
    open_about(window);
    CHECK(child<QTabWidget>(window, u"about_tabs").count() == 4);
    // packaging/licenses/THIRD-PARTY.md, LICENSE, packaging/licenses/LGPL-3.0.txt
    auto const third_party = child<QTextBrowser>(window, u"third_party_text").toPlainText();
    CHECK(third_party.contains(u"Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. (2010)"_s));
    CHECK(third_party.contains(u"https://download.qt.io/archive/qt/"_s));
    CHECK(child<QTextBrowser>(window, u"gpl_text").toPlainText().contains(u"GNU GENERAL PUBLIC LICENSE"_s));
    auto const lgpl = child<QTextBrowser>(window, u"lgpl_text").toPlainText();
    CHECK(lgpl.contains(u"GNU LESSER GENERAL PUBLIC LICENSE"_s));
    CHECK(lgpl.contains(u"Version 3, 29 June 2007"_s));
}
