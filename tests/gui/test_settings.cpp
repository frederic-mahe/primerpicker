#include "gui/colour_scheme.hpp"
#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QCheckBox>
#include <QLabel>
#include <QLineEdit>
#include <QSettings>
#include <QString>
#include <QStringView>
#include <QTableView>
#include <QTemporaryDir>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::ColourScheme;
using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

void direct_entry(MainWindow& window, QString const& primer) {
    child<QAction>(window, u"direct_entry").trigger();
    auto& field = child<QLineEdit>(window, u"direct_entry_field");
    QTest::keyClicks(&field, primer);
    QTest::keyClick(&field, Qt::Key_Return);
}

auto settings_column(MainWindow const& window) -> QList<QString> {
    auto const& model = *child<QTableView>(window, u"criteria_table").model();
    constexpr auto setting_column = 1;
    auto settings = QList<QString>{};
    for (auto row = 0; row < model.rowCount(); ++row) {
        settings.append(model.data(model.index(row, setting_column)).toString());
    }
    return settings;
}

}  // namespace

// Reference: PRIMER.exe, One Primer report of the M13 primer with
// config/PD2_REPORT_NARROW.CFG as PD2.CFG
// (tests/legacy/reports/inspect/GTAAAACGACGGCCAGT_NARROW.PRN, recorded
// 2026-10-20)
TEST_CASE("a pd2.cfg file replaces the built-in criteria", "[gui][settings]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_settings(legacy_file(u"config/PD2_REPORT_NARROW.CFG")).has_value());
    CHECK(child<QLabel>(window, u"settings_label").text() == u"Settings: PD2_REPORT_NARROW.CFG"_s);
    direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    CHECK(settings_column(window) ==
          QList<QString>{u"Min 40, Max 65"_s, u"Min 50, Max 70"_s, u"Energy cutoff -3.5 kcal"_s,
                         u"Reject >= 4 matches at 3' end"_s, u"Reject >= 6 adj homol bases"_s,
                         u"Reject >= 4 base runs"_s, u"Reject >= 2 G or C at 3' end"_s});
}

// Reference: config/PD2_LENGTH25.CFG holds the built-in defaults with a
// primer length of 25 (tests/legacy/README.md, "Settings file
// pd2.cfg"); From Molecule makes primers of that length
TEST_CASE("From Molecule uses the primer length of the settings file", "[gui][settings]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_settings(legacy_file(u"config/PD2_LENGTH25.CFG")).has_value());
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    child<QAction>(window, u"from_molecule").trigger();
    auto& field = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&field, u"51"_s);
    QTest::keyClick(&field, Qt::Key_Return);
    CHECK(child<QLabel>(window, u"length_label").text() == u"25"_s);
}

TEST_CASE("a file that is not a pd2.cfg is refused, and the criteria are kept", "[gui][settings]") {
    auto window = MainWindow{};
    CHECK(window.load_settings(legacy_file(u"pairs/TEST600.GB")).has_value());
    CHECK(child<QLabel>(window, u"settings_label").text().isEmpty());
    direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    CHECK(settings_column(window).at(0) == u"Min 50, Max 60"_s);
}

TEST_CASE("colour schemes and their names in the settings", "[gui][colour scheme]") {
    CHECK(primerpicker::gui::colour_scheme_name(ColourScheme::system) == u"system"_s);
    CHECK(primerpicker::gui::colour_scheme_name(ColourScheme::light) == u"light"_s);
    CHECK(primerpicker::gui::colour_scheme_name(ColourScheme::dark) == u"dark"_s);
    CHECK(primerpicker::gui::colour_scheme_from_name(u"dark"_s) == ColourScheme::dark);
    CHECK(primerpicker::gui::colour_scheme_from_name(u"light"_s) == ColourScheme::light);
    // anything else: follow the system
    CHECK(primerpicker::gui::colour_scheme_from_name(u"blue"_s) == ColourScheme::system);
    CHECK(primerpicker::gui::colour_scheme_from_name({}) == ColourScheme::system);
}

TEST_CASE("the colour scheme follows the system by default", "[gui][colour scheme]") {
    auto const window = MainWindow{};
    CHECK(child<QAction>(window, u"colour_scheme_system").isChecked());
    CHECK_FALSE(child<QAction>(window, u"colour_scheme_dark").isChecked());
}

TEST_CASE("the chosen colour scheme is kept in the settings", "[gui][colour scheme]") {
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto settings = QSettings{directory.filePath(u"primerpicker.ini"_s), QSettings::IniFormat};
    {
        auto window = MainWindow{settings};
        child<QAction>(window, u"colour_scheme_dark").trigger();
        CHECK(child<QAction>(window, u"colour_scheme_dark").isChecked());
        CHECK_FALSE(child<QAction>(window, u"colour_scheme_system").isChecked());
    }
    CHECK(settings.value(u"appearance/colour_scheme"_s).toString() == u"dark"_s);
    auto const window = MainWindow{settings};
    CHECK(child<QAction>(window, u"colour_scheme_dark").isChecked());
    primerpicker::gui::apply_colour_scheme(ColourScheme::system);  // leave the application as found
}
