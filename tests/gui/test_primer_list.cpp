#include "gui/main_window.hpp"
#include "gui_primer_helpers.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QCheckBox>
#include <QDockWidget>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QString>
#include <QStringView>
#include <QTableView>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;
using primerpicker::gui::test::direct_entry;
using primerpicker::gui::test::save_direct_entry;
using primerpicker::gui::test::table_rows;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

}  // namespace

TEST_CASE("Save to List adds the active primer to the Primer List", "[gui][primer list]") {
    auto window = MainWindow{};
    CHECK(table_rows(window, u"primer_list").isEmpty());
    save_direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    save_direct_entry(window, u"AACAGCTATGACCATG"_s);
    CHECK(table_rows(window, u"primer_list") ==
          QList<QList<QString>>{{u"1"_s, u"GTAAAACGACGGCCAGT"_s, u"17"_s, u"Direct Entry"_s},
                                {u"2"_s, u"AACAGCTATGACCATG"_s, u"16"_s, u"Direct Entry"_s}});
}

TEST_CASE("the Primer List shows where a primer comes from", "[gui][primer list]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    for (auto const complement : {false, true}) {
        child<QAction>(window, u"from_molecule").trigger();
        child<QCheckBox>(window, u"complement_box").setChecked(complement);
        auto& position = child<QLineEdit>(window, u"position_field");
        QTest::keyClicks(&position, complement ? u"400"_s : u"51"_s);
        QTest::keyClick(&position, Qt::Key_Return);
        child<QAction>(window, u"save_to_list").trigger();
    }
    REQUIRE_FALSE(window.load_primer(legacy_file(u"primers/M13EXT.PD2")).has_value());
    child<QAction>(window, u"save_to_list").trigger();
    // TEST600.GB at 51 and at 400 (complement strand): pairs/pairs.tsv
    CHECK(table_rows(window, u"primer_list") ==
          QList<QList<QString>>{{u"1"_s, u"CTCAAAGGCTCAGACATCAG"_s, u"20"_s, u"TEST600 at 51"_s},
                                {u"2"_s, u"TGTATCCAGGTAATGCAAGG"_s, u"20"_s, u"TEST600 at 400 (complement)"_s},
                                {u"3"_s, u"GAATTC-GTAAAACGACGGCCAGT"_s, u"23"_s, u"M13EXT.PD2"_s}});
}

TEST_CASE("Save to List without a primer asks for one", "[gui][primer list]") {
    auto window = MainWindow{};
    child<QAction>(window, u"save_to_list").trigger();
    CHECK(table_rows(window, u"primer_list").isEmpty());
    CHECK(child<QLabel>(window, u"no_primer_label").isVisibleTo(&window));
}

TEST_CASE("a primer of the list becomes the active primer", "[gui][primer list]") {
    auto window = MainWindow{};
    save_direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    save_direct_entry(window, u"AACAGCTATGACCATG"_s);
    auto& list = child<QTableView>(window, u"primer_list");
    list.setCurrentIndex(list.model()->index(0, 0));
    QTest::keyClick(&list, Qt::Key_Return);
    CHECK(child<QLabel>(window, u"sequence_label").text() == u"5'-GTAAAACGACGGCCAGT-3'"_s);
}

TEST_CASE("View -> Primer List shows or hides the list", "[gui][primer list]") {
    auto window = MainWindow{};
    window.show();  // the menu item follows the dock's visibility
    auto& toggle = child<QAction>(window, u"view_primer_list");
    auto const& dock = child<QDockWidget>(window, u"primer_list_dock");
    CHECK(dock.isVisibleTo(&window));
    toggle.trigger();
    CHECK_FALSE(dock.isVisibleTo(&window));
    toggle.trigger();
    CHECK(dock.isVisibleTo(&window));
}

TEST_CASE("Direct Entry alone does not add to the list", "[gui][primer list]") {
    auto window = MainWindow{};
    direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    CHECK(table_rows(window, u"primer_list").isEmpty());
}
