#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QCheckBox>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QString>
#include <QStringView>
#include <QTableView>
#include <QTest>
#include <QWidget>

#include <catch2/catch_test_macros.hpp>

#include <stdexcept>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

enum struct Complement : bool { no, yes };

// From Molecule at POSITION of TEST600, saved with S (entry 1 of the
// Primer List, active)
void saved_primer(MainWindow& window, QString const& position, Complement const complement = Complement::no) {
    if (window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value()) {
        throw std::runtime_error{"cannot load TEST600.GB"};
    }
    child<QAction>(window, u"from_molecule").trigger();
    child<QCheckBox>(window, u"complement_box").setChecked(complement == Complement::yes);
    auto& field = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&field, position);
    QTest::keyClick(&field, Qt::Key_Return);
    QTest::keyClick(&child<QWidget>(window, u"picking_view"), Qt::Key_S);
}

auto list_rows(MainWindow const& window) -> QList<QList<QString>> {
    auto const& model = *child<QTableView>(window, u"primer_list").model();
    auto rows = QList<QList<QString>>{};
    for (auto row = 0; row < model.rowCount(); ++row) {
        auto cells = QList<QString>{};
        for (auto column = 0; column < model.columnCount(); ++column) {
            cells.append(model.data(model.index(row, column)).toString());
        }
        rows.append(cells);
    }
    return rows;
}

}  // namespace

// Reference: PRIMER.exe, Modify -> Length on the primer at 51 of TEST600,
// saved as entry 1 (tests/legacy/modify_length/, recorded 2026-10-21):
// 25 changes entry 1 in place to CTCAAAGGCTCAGACATCAGGGCTC
TEST_CASE("Modify -> Length changes the active primer and its Primer List entry", "[gui][length]") {
    auto window = MainWindow{};
    saved_primer(window, u"51"_s);
    REQUIRE_FALSE(window.change_length(25).has_value());
    CHECK(child<QLabel>(window, u"sequence_label").text() == u"5'-CTCAAAGGCTCAGACATCAGGGCTC-3'"_s);
    CHECK(list_rows(window) ==
          QList<QList<QString>>{{u"1"_s, u"CTCAAAGGCTCAGACATCAGGGCTC"_s, u"25"_s, u"TEST600 at 51"_s}});
}

// Same reference: 400(C), 25 bases: TGTATCCAGGTAATGCAAGGATACC
TEST_CASE("Modify -> Length on the complement strand keeps the 5' end", "[gui][length]") {
    auto window = MainWindow{};
    saved_primer(window, u"400"_s, Complement::yes);
    REQUIRE_FALSE(window.change_length(25).has_value());
    CHECK(child<QLabel>(window, u"sequence_label").text() == u"5'-TGTATCCAGGTAATGCAAGGATACC-3'"_s);
}

// Same reference: 0 and 30 bases at 575 are refused with PRIMER.exe's
// messages; the primer is kept
TEST_CASE("Modify -> Length refuses no base and a primer off the molecule", "[gui][length]") {
    auto window = MainWindow{};
    saved_primer(window, u"575"_s);
    CHECK(window.change_length(0) == u"Primer must have length!"_s);
    CHECK(window.change_length(30) == u"Primer Too Long - Extends outside molecule"_s);
    CHECK(child<QLabel>(window, u"length_label").text() == u"20"_s);
}

TEST_CASE("Modify -> Length needs a primer made from the loaded molecule", "[gui][length]") {
    auto window = MainWindow{};
    CHECK(window.change_length(25).has_value());
    REQUIRE_FALSE(window.load_primer(legacy_file(u"primers/M13EXT.PD2")).has_value());
    CHECK(window.change_length(25).has_value());
}
