#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QDialog>
#include <QFile>
#include <QIODevice>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QString>
#include <QStringView>
#include <QTableView>
#include <QTemporaryDir>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

// The field of a sequence dialog, and the menu item that opens it
struct SequenceInput {
    QStringView action;
    QStringView field;
};

constexpr auto direct_entry = SequenceInput{.action = u"direct_entry", .field = u"direct_entry_field"};
constexpr auto extension = SequenceInput{.action = u"modify_extension", .field = u"extension_field"};

// types SEQUENCE into the field, then Enter
void type_in(MainWindow& window, SequenceInput const input, QString const& sequence) {
    child<QAction>(window, input.action).trigger();
    auto& line = child<QLineEdit>(window, input.field);
    QTest::keyClicks(&line, sequence);
    QTest::keyClick(&line, Qt::Key_Return);
}

auto label_text(MainWindow const& window, QStringView const name) -> QString {
    return child<QLabel>(window, name).text();
}

auto column(MainWindow const& window, int const index) -> QList<QString> {
    auto const& model = *child<QTableView>(window, u"criteria_table").model();
    auto values = QList<QString>{};
    for (auto row = 0; row < model.rowCount(); ++row) {
        values.append(model.data(model.index(row, index)).toString());
    }
    return values;
}

auto read_all(QString const& path) -> QByteArray {
    auto file = QFile{path};
    REQUIRE(file.open(QIODevice::ReadOnly));
    return file.readAll();
}

// the text after the first line (the first line holds the report date)
auto after_first_line(QByteArray const& report) -> QByteArray {
    return report.mid(report.indexOf('\n') + 1);
}

constexpr auto actual_column = 2;
constexpr auto meets_column = 3;

}  // namespace

// Reference: PRIMER.exe, One Primer report of the M13 primer with the
// 5' extension CCCC (tests/legacy/reports/inspect/
// CCCC_GTAAAACGACGGCCAGT.PRN, recorded 2026-10-20), and the hairpin of
// the Inspect screens (tests/legacy/README.md, "5' extension and the
// criteria", recorded 2026-10-18: -2.7)
TEST_CASE("Modify -> 5'-extension adds an extension to the active primer", "[gui][extension]") {
    auto window = MainWindow{};
    type_in(window, direct_entry, u"GTAAAACGACGGCCAGT"_s);
    type_in(window, extension, u"CCCC"_s);
    CHECK_FALSE(child<QDialog>(window, u"extension_dialog").isVisible());
    CHECK(label_text(window, u"sequence_label") == u"5'-CCCC-GTAAAACGACGGCCAGT-3'"_s);
    CHECK(label_text(window, u"length_label") == u"21"_s);
    CHECK(label_text(window, u"extension_label") == u"CCCC"_s);
    CHECK(label_text(window, u"composition_label") == u"A: 6  C: 8  G: 5  T: 2"_s);
    CHECK(column(window, actual_column) ==
          QList<QString>{u"61"_s, u"64"_s, u"--"_s, u"2"_s, u"4"_s, u"4"_s, u"0"_s});
    CHECK(column(window, meets_column) ==
          QList<QString>{u"NO"_s, u"YES"_s, u"NO"_s, u"YES"_s, u"YES"_s, u"NO"_s, u"YES"_s});
    CHECK(label_text(window, u"hairpin_label") == u"-2.7 kcal"_s);
}

TEST_CASE("the extension field opens with the current extension, and an empty field removes it",
          "[gui][extension]") {
    auto window = MainWindow{};
    type_in(window, direct_entry, u"GTAAAACGACGGCCAGT"_s);
    type_in(window, extension, u"CCCC"_s);
    child<QAction>(window, u"modify_extension").trigger();
    auto& field = child<QLineEdit>(window, u"extension_field");
    CHECK(field.text() == u"CCCC"_s);
    field.clear();
    QTest::keyClick(&field, Qt::Key_Return);
    CHECK(label_text(window, u"sequence_label") == u"5'-GTAAAACGACGGCCAGT-3'"_s);
    CHECK(label_text(window, u"extension_label") == u"--none--"_s);
}

TEST_CASE("Modify -> 5'-extension without a primer asks for one", "[gui][extension]") {
    auto window = MainWindow{};
    child<QAction>(window, u"modify_extension").trigger();
    CHECK_FALSE(child<QDialog>(window, u"extension_dialog").isVisible());
    CHECK(child<QLabel>(window, u"no_primer_label").isVisibleTo(&window));
}

// Reference: primers/M13EXT.PD2, written by PRIMER.exe's Save Active
// Primer (tests/legacy/README.md, "Primer files .PD2", recorded
// 2026-10-18): Direct Entry, 5' extension GAATTC, primer
// GTAAAACGACGGCCAGT, built-in criteria
TEST_CASE("Primer/Create -> Retrieve from disk loads a .PD2 primer", "[gui][pd2]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_primer(legacy_file(u"primers/M13EXT.PD2")).has_value());
    CHECK(label_text(window, u"sequence_label") == u"5'-GAATTC-GTAAAACGACGGCCAGT-3'"_s);
    CHECK(label_text(window, u"extension_label") == u"GAATTC"_s);
    CHECK(column(window, 1).at(0) == u"Min 50, Max 60"_s);
}

// Same reference: primers/COMP.PD2, from TEST600.GB at 400, complement
// strand: TGTATCCAGGTAATGCAAGG
TEST_CASE("a .PD2 primer made from a molecule", "[gui][pd2]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_primer(legacy_file(u"primers/COMP.PD2")).has_value());
    CHECK(label_text(window, u"sequence_label") == u"5'-TGTATCCAGGTAATGCAAGG-3'"_s);
}

TEST_CASE("a file that is not a .PD2 primer is refused, and the active primer kept", "[gui][pd2]") {
    auto window = MainWindow{};
    type_in(window, direct_entry, u"GTAAAACGACGGCCAGT"_s);
    CHECK(window.load_primer(legacy_file(u"pairs/TEST600.GB")).has_value());
    CHECK(label_text(window, u"sequence_label") == u"5'-GTAAAACGACGGCCAGT-3'"_s);
}

// Reference: PRIMER.exe, Print to Disk of the Primer Analysis screen of
// the M13 primer (tests/legacy/reports/INSPECT.PRN, recorded
// 2026-10-19), byte for byte but for the date line
TEST_CASE("Print to Disk writes PRIMER.exe's One Primer report", "[gui][report]") {
    auto window = MainWindow{};
    type_in(window, direct_entry, u"GTAAAACGACGGCCAGT"_s);
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"INSPECT.PRN"_s);
    REQUIRE_FALSE(window.print_report_to(path).has_value());
    auto const report = read_all(path);
    CHECK(report.startsWith("Primer Design Report:  One Primer"_ba));
    CHECK(after_first_line(report) == after_first_line(read_all(legacy_file(u"reports/INSPECT.PRN"))));
}

TEST_CASE("Print to Disk without a primer writes nothing", "[gui][report]") {
    auto window = MainWindow{};
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"INSPECT.PRN"_s);
    CHECK(window.print_report_to(path).has_value());
    CHECK_FALSE(QFile::exists(path));
}
