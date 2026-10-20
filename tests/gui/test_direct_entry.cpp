#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QDialog>
#include <QDialogButtonBox>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QPushButton>
#include <QString>
#include <QStringView>
#include <QAbstractItemModel>
#include <QTableView>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto open_direct_entry(MainWindow& window) -> QLineEdit& {
    child<QAction>(window, u"direct_entry").trigger();
    return child<QLineEdit>(window, u"direct_entry_field");
}

void direct_entry(MainWindow& window, QString const& primer) {
    auto& field = open_direct_entry(window);
    QTest::keyClicks(&field, primer);
    QTest::keyClick(&field, Qt::Key_Return);
}

auto row_texts(QTableView const& table, int const row) -> QList<QString> {
    auto const& model = *table.model();
    auto texts = QList<QString>{};
    for (auto column = 0; column < model.columnCount(); ++column) {
        texts.append(model.data(model.index(row, column)).toString());
    }
    return texts;
}

auto label_text(MainWindow const& window, QStringView const name) -> QString {
    return child<QLabel>(window, name).text();
}

}  // namespace

TEST_CASE("Primer/Create -> Direct Entry opens a dialog with an empty field", "[gui][direct entry]") {
    auto window = MainWindow{};
    auto const& field = open_direct_entry(window);
    CHECK(child<QDialog>(window, u"direct_entry_dialog").isVisible());
    CHECK(field.text().isEmpty());
}

// Reference: PRIMER.exe's Direct Entry field (tests/legacy/README.md,
// "Primer Tm", recorded 2026-10-17): it holds 50 bases and does not
// insert letters other than A, C, G and T
TEST_CASE("the Direct Entry field inserts A, C, G and T only", "[gui][direct entry]") {
    auto window = MainWindow{};
    auto& field = open_direct_entry(window);
    QTest::keyClicks(&field, u"ACGTNRUacgt"_s);
    CHECK(field.text() == u"ACGTACGT"_s);
}

TEST_CASE("the Direct Entry field holds 50 bases", "[gui][direct entry]") {
    auto window = MainWindow{};
    auto& field = open_direct_entry(window);
    constexpr auto typed_length = 60;
    QTest::keyClicks(&field, QString(typed_length, u'A'));
    CHECK(field.text() == QString(50, u'A'));
}

TEST_CASE("Direct Entry needs two bases", "[gui][direct entry]") {
    auto window = MainWindow{};
    auto& field = open_direct_entry(window);
    auto const& ok_button = *child<QDialogButtonBox>(window, u"direct_entry_buttons").button(QDialogButtonBox::Ok);
    QTest::keyClicks(&field, u"A"_s);
    CHECK_FALSE(ok_button.isEnabled());
    QTest::keyClick(&field, Qt::Key_Return);
    CHECK(child<QDialog>(window, u"direct_entry_dialog").isVisible());
    QTest::keyClicks(&field, u"C"_s);
    CHECK(ok_button.isEnabled());
}

TEST_CASE("Escape leaves Direct Entry without a primer", "[gui][direct entry]") {
    auto window = MainWindow{};
    auto& field = open_direct_entry(window);
    QTest::keyClicks(&field, u"GTAAAACGACGGCCAGT"_s);
    QTest::keyClick(&field, Qt::Key_Escape);
    CHECK_FALSE(child<QDialog>(window, u"direct_entry_dialog").isVisible());
    CHECK(child<QLabel>(window, u"no_primer_label").isVisibleTo(&window));
}

// Reference: PRIMER.exe's message when there is no active primer (from
// the strings of examples/PRIMER.exe)
TEST_CASE("Analyze -> Active Primer without a primer asks for one", "[gui][analysis]") {
    auto window = MainWindow{};
    child<QAction>(window, u"analyze_active_primer").trigger();
    auto const& label = child<QLabel>(window, u"no_primer_label");
    CHECK(label.isVisibleTo(&window));
    CHECK(label.text().startsWith(u"Please create or retrieve a primer:"_s));
}

// Reference: PRIMER.exe, One Primer report of the M13 primer typed with
// Direct Entry, built-in defaults (tests/legacy/reports/INSPECT.PRN,
// recorded 2026-10-19)
TEST_CASE("Direct Entry shows the Primer Analysis of the typed primer", "[gui][analysis]") {
    auto window = MainWindow{};
    direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    CHECK_FALSE(child<QDialog>(window, u"direct_entry_dialog").isVisible());
    CHECK(child<QWidget>(window, u"analysis_view").isVisibleTo(&window));
    CHECK(label_text(window, u"sequence_label") == u"5'-GTAAAACGACGGCCAGT-3'"_s);
    CHECK(label_text(window, u"length_label") == u"17"_s);
    CHECK(label_text(window, u"extension_label") == u"--none--"_s);
    CHECK(label_text(window, u"composition_label") == u"A: 6  C: 4  G: 5  T: 2"_s);
    auto const& table = child<QTableView>(window, u"criteria_table");
    REQUIRE(table.model()->rowCount() == 7);
    CHECK(row_texts(table, 0) == QList<QString>{u"% GC"_s, u"Min 50, Max 60"_s, u"52"_s, u"YES"_s});
    CHECK(row_texts(table, 1) == QList<QString>{u"Tm °C"_s, u"Min 55, Max 80"_s, u"64"_s, u"YES"_s});
    CHECK(row_texts(table, 2) == QList<QString>{u"No Hairpins"_s, u"Energy cutoff 0.0 kcal"_s, u"--"_s, u"YES"_s});
    CHECK(row_texts(table, 3) ==
          QList<QString>{u"No 3' Dimers"_s, u"Reject >= 3 matches at 3' end"_s, u"2"_s, u"YES"_s});
    CHECK(row_texts(table, 4) == QList<QString>{u"No Dimers"_s, u"Reject >= 7 adj homol bases"_s, u"4"_s, u"YES"_s});
    CHECK(row_texts(table, 5) == QList<QString>{u"No Runs"_s, u"Reject >= 3 base runs"_s, u"4"_s, u"NO"_s});
    CHECK(row_texts(table, 6) ==
          QList<QString>{u"No 3'GC runs"_s, u"Reject >= 3 G or C at 3' end"_s, u"0"_s, u"YES"_s});
    CHECK(label_text(window, u"verdict_label") == u"NO"_s);
}

// Reference: PRIMER.exe, Direct Entry -> Inspect and H (Hairpins
// screen), built-in defaults (tests/legacy/inspect/inspect.tsv and
// tm/direct_entry_tm.tsv, recorded 2026-10-17): %GC 66, Tm 62, 3'
// dimers 3, dimers 3, runs 2, 3' GC runs 2, hairpin -.4 kcal
TEST_CASE("the Primer Analysis view of a primer with a hairpin", "[gui][analysis]") {
    auto window = MainWindow{};
    direct_entry(window, u"CGTCCAGCACCTTGG"_s);
    auto const& table = child<QTableView>(window, u"criteria_table");
    auto actual = QList<QString>{};
    for (auto row = 0; row < table.model()->rowCount(); ++row) {
        actual.append(row_texts(table, row).at(2));
    }
    CHECK(actual == QList<QString>{u"66"_s, u"62"_s, u"--"_s, u"3"_s, u"3"_s, u"2"_s, u"2"_s});
    CHECK(row_texts(table, 2).at(3) == u"NO"_s);
    CHECK(label_text(window, u"hairpin_label") == u"-0.4 kcal"_s);
    CHECK(label_text(window, u"verdict_label") == u"NO"_s);
}

TEST_CASE("Analyze -> Active Primer shows the active primer again", "[gui][analysis]") {
    auto window = MainWindow{};
    direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    child<QAction>(window, u"analyze_active_primer").trigger();
    CHECK(child<QWidget>(window, u"analysis_view").isVisibleTo(&window));
    CHECK(label_text(window, u"sequence_label") == u"5'-GTAAAACGACGGCCAGT-3'"_s);
}
