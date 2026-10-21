#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QApplication>
#include <QByteArray>
#include <QCheckBox>
#include <QDialog>
#include <QDir>
#include <QFile>
#include <QIODevice>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QListWidget>
#include <QRegularExpression>
#include <QString>
#include <QStringView>
#include <QTableView>
#include <QTemporaryDir>
#include <QTest>
#include <QWidget>
#include <QtEnvironmentVariables>

#include <catch2/catch_test_macros.hpp>

#include <stdexcept>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

// From Molecule at POSITION of TEST600 (top strand), saved with S, then
// Modify -> Mutagenize
void mutagenize_test600(MainWindow& window, QString const& position) {
    if (window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value()) {
        throw std::runtime_error{"cannot load TEST600.GB"};
    }
    child<QAction>(window, u"from_molecule").trigger();
    child<QCheckBox>(window, u"complement_box").setChecked(false);
    auto& field = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&field, position);
    QTest::keyClick(&field, Qt::Key_Return);
    QTest::keyClick(&child<QWidget>(window, u"picking_view"), Qt::Key_S);
    REQUIRE_FALSE(window.mutagenize().has_value());
}

auto view(MainWindow& window) -> QWidget& {
    return child<QWidget>(window, u"mutagenize_view");
}

void keys(MainWindow& window, QString const& typed) {
    QTest::keyClicks(&view(window), typed);
}

void right(MainWindow& window, int const count) {
    for (auto index = 0; index < count; ++index) {
        QTest::keyClick(&view(window), Qt::Key_Right);
    }
}

auto original_line(MainWindow const& window) -> QString {
    return child<QLabel>(window, u"mutagenize_original").text();
}

auto mutated_line(MainWindow const& window) -> QString {
    return child<QLabel>(window, u"mutagenize_mutated").text();
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

// a report without its dates (one per page header)
auto body(QString const& path) -> QByteArray {
    auto file = QFile{path};
    REQUIRE(file.open(QIODevice::ReadOnly));
    auto text = file.readAll();
    static auto const date = QRegularExpression{u"Report Date: [0-9/]+"_s};
    return QString::fromLatin1(text).replace(date, u"Report Date: DATE"_s).toLatin1();
}

}  // namespace

// Reference: PRIMER.exe, Mutagenize on TEST600 at 52, no edit
// (tests/legacy/mutagenize/mutagenize.tsv, recorded 2026-10-19)
TEST_CASE("Modify -> Mutagenize shows the original and mutated values", "[gui][mutagenize]") {
    auto window = MainWindow{};
    mutagenize_test600(window, u"52"_s);
    CHECK(view(window).isVisibleTo(&window));
    CHECK(original_line(window) == u"ORIGINAL: %GC: 50    Tm°C: 69    Hairpin: OK  Dimers: OK  Runs: !"_s);
    CHECK(mutated_line(window) == u"MUTATED:  %GC: 50    Tm°C: 69    Hairpin: OK  Dimers: OK  Runs: !"_s);
}

// Reference: PRIMER.exe, Mutagenize on TEST600 at 51 with the keys
// "a right right right t right g" (mutagenize.tsv): typed bases replace
// the highlighted one, which does not move
TEST_CASE("typed bases replace the highlighted base", "[gui][mutagenize]") {
    auto window = MainWindow{};
    mutagenize_test600(window, u"51"_s);
    keys(window, u"a"_s);
    right(window, 3);
    keys(window, u"t"_s);
    right(window, 1);
    keys(window, u"g"_s);
    CHECK(original_line(window) == u"ORIGINAL: %GC: 50    Tm°C: 66    Hairpin: OK  Dimers: OK  Runs: !"_s);
    CHECK(mutated_line(window) == u"MUTATED:  %GC: 50    Tm°C: 65    Hairpin: !   Dimers: !   Runs: OK"_s);
}

// Reference: PRIMER.exe (tests/legacy/mutagenize_enter/,
// mutagenize_f1/, recorded 2026-10-20 and 21): S saves the mutated
// primer as a new entry of the Primer List ("Primer Saved"), which is
// then active; Mutagenize on it: "Primer is already Mutated"; Modify ->
// Position: "Mutated primers cannot be moved"
TEST_CASE("S saves the mutated primer to the Primer List", "[gui][mutagenize]") {
    auto window = MainWindow{};
    mutagenize_test600(window, u"51"_s);
    keys(window, u"g"_s);
    right(window, 1);
    keys(window, u"a"_s);
    right(window, 1);
    keys(window, u"a"_s);
    right(window, 1);
    keys(window, u"t"_s);
    right(window, 1);
    keys(window, u"t"_s);
    right(window, 1);
    keys(window, u"c"_s);
    QTest::keyClick(&view(window), Qt::Key_S);
    CHECK(child<QLabel>(window, u"mutagenize_message").text() == u"Primer Saved"_s);
    CHECK(list_rows(window) ==
          QList<QList<QString>>{{u"1"_s, u"CTCAAAGGCTCAGACATCAG"_s, u"20"_s, u"TEST600 at 51"_s},
                                {u"2"_s, u"GAATTCGGCTCAGACATCAG"_s, u"20"_s, u"TEST600 at 51 (mutated)"_s}});
    QTest::keyClick(&view(window), Qt::Key_Escape);
    CHECK(window.mutagenize() == u"Primer is already Mutated"_s);
    CHECK(window.modify_position() == u"Mutated primers cannot be moved"_s);
}

// Reference: PRIMER.exe, GAATTC typed over the first six bases, then I
// and Print to Disk (tests/legacy/mutagenize_enter/
// TEST600_51_GAATTC_i_report.PRN, recorded 2026-10-21)
TEST_CASE("I, then Print to Disk: the One Primer report of the mutated primer", "[gui][mutagenize][report]") {
    auto window = MainWindow{};
    mutagenize_test600(window, u"51"_s);
    for (auto const* const base : {u"g", u"a", u"a", u"t", u"t"}) {
        keys(window, QString::fromUtf16(base));
        right(window, 1);
    }
    keys(window, u"c"_s);
    QTest::keyClick(&view(window), Qt::Key_I);
    CHECK(child<QWidget>(window, u"analysis_view").isVisibleTo(&window));
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"REPORT"_s);
    REQUIRE_FALSE(window.print_report_to(path).has_value());
    CHECK(body(path) == body(legacy_file(u"mutagenize_enter/TEST600_51_GAATTC_i_report.PRN")));
}

// Reference: PRIMER.exe opens its enzyme list with F1 and the codon
// table with F2 (mutagenize_f1/, mutagenize/codon_table.tsv); here F1
// lists PrimerPicker's enzymes
TEST_CASE("F1 lists the enzymes, F2 the codon table", "[gui][mutagenize]") {
    auto window = MainWindow{};
    mutagenize_test600(window, u"51"_s);
    QTest::keyClick(&view(window), Qt::Key_F1);
    auto const& enzymes = child<QDialog>(window, u"mutagenize_enzymes");
    CHECK(enzymes.isVisible());
    CHECK(child<QListWidget>(window, u"mutagenize_enzymes_list").count() > 500);
    child<QDialog>(window, u"mutagenize_enzymes").reject();
    QTest::keyClick(&view(window), Qt::Key_F2);
    auto const& codons = child<QListWidget>(window, u"mutagenize_codons_list");
    CHECK(codons.count() == 20);
    // CTC, the first codon: Leu
    CHECK(codons.currentItem()->text().startsWith(u"Leu"_s));
}

TEST_CASE("Mutagenize needs a primer made from the loaded molecule", "[gui][mutagenize]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_primer(legacy_file(u"primers/M13EXT.PD2")).has_value());
    CHECK(window.mutagenize().has_value());
}

// Hidden test (not run by CTest), as tests/gui/test_screenshots.cpp:
//   PRIMERPICKER_SCREENSHOT_DIR=DIR build/debug/tests/gui/gui_tests '[.mutagenize_screenshots]'
TEST_CASE("a screenshot of the Mutagenize view", "[.mutagenize_screenshots]") {
    auto const directory = qEnvironmentVariable("PRIMERPICKER_SCREENSHOT_DIR");
    REQUIRE_FALSE(directory.isEmpty());
    QApplication::setStyle(u"Fusion"_s);
    auto window = MainWindow{};
    constexpr auto width = 900;
    constexpr auto height = 560;
    window.resize(width, height);
    window.show();
    mutagenize_test600(window, u"51"_s);
    keys(window, u"g"_s);
    right(window, 3);
    QApplication::processEvents();
    CHECK(window.grab().save(QDir{directory}.filePath(u"mutagenize_light.png"_s)));
}
