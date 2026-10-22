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

// From Molecule at POSITION of TEST600 on the complement strand, saved
// with S, then Modify -> Mutagenize
void mutagenize_complement(MainWindow& window, QString const& position) {
    if (window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value()) {
        throw std::runtime_error{"cannot load TEST600.GB"};
    }
    child<QAction>(window, u"from_molecule").trigger();
    child<QCheckBox>(window, u"complement_box").setChecked(true);
    auto& field = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&field, position);
    QTest::keyClick(&field, Qt::Key_Return);
    QTest::keyClick(&child<QWidget>(window, u"picking_view"), Qt::Key_S);
    REQUIRE_FALSE(window.mutagenize().has_value());
}

auto view(MainWindow& window) -> QWidget& {
    return child<QWidget>(window, u"mutagenize_view");
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

// Reference: PRIMER.exe, Mutagenize on TEST600 at 400(C), no edit
// (tests/legacy/complement/mut_c400.screens.txt, recorded 2026-10-22):
// the primer's own strand, "COMPLEMENTARY" in the header
TEST_CASE("Mutagenize shows a complement-strand primer on its own strand", "[gui][mutagenize]") {
    auto window = MainWindow{};
    mutagenize_complement(window, u"400"_s);
    CHECK(view(window).isVisibleTo(&window));
    CHECK(child<QLabel>(window, u"mutagenize_header").text().contains(u"COMPLEMENTARY"_s));
    auto const screen = child<QLabel>(window, u"mutagenize_screen").text();
    CHECK(screen.contains(u"TCATCTGGCGTAAGTGTCAG"_s));  // the molecule line's first bases
    CHECK(child<QLabel>(window, u"mutagenize_original").text() ==
          u"ORIGINAL: %GC: 45!   Tm°C: 66    Hairpin: !   Dimers: OK  Runs: OK"_s);
}

// Reference: PRIMER.exe, a over base 3 and g over base 7 of 400(C), then
// S (complement/mut_c400_edit.screens.txt): MUTATED Tm 69, Hairpin OK,
// "Primer Saved"
TEST_CASE("Mutagenize saves a mutated complement-strand primer", "[gui][mutagenize]") {
    auto window = MainWindow{};
    mutagenize_complement(window, u"400"_s);
    QTest::keyClick(&view(window), Qt::Key_Right);
    QTest::keyClick(&view(window), Qt::Key_Right);
    QTest::keyClicks(&view(window), u"a"_s);
    for (auto index = 0; index < 4; ++index) {
        QTest::keyClick(&view(window), Qt::Key_Right);
    }
    QTest::keyClicks(&view(window), u"g"_s);
    CHECK(child<QLabel>(window, u"mutagenize_mutated").text() ==
          u"MUTATED:  %GC: 45!   Tm°C: 69    Hairpin: OK  Dimers: OK  Runs: OK"_s);
    QTest::keyClick(&view(window), Qt::Key_S);
    CHECK(child<QLabel>(window, u"mutagenize_message").text() == u"Primer Saved"_s);
    auto const rows = list_rows(window);
    REQUIRE(rows.size() == 2);
    CHECK(rows.at(1).at(1) == u"TGAATCGAGGTAATGCAAGG"_s);
}
