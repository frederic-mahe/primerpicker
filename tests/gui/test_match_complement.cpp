#include "gui/main_window.hpp"

#include "gui_primer_helpers.hpp"
#include "gui_test_helpers.hpp"

#include "core/sequence/strand.hpp"

#include <QAction>
#include <QByteArray>
#include <QCheckBox>
#include <QDialog>
#include <QFile>
#include <QIODevice>
#include <QLineEdit>
#include <QList>
#include <QRegularExpression>
#include <QString>
#include <QStringView>
#include <QTemporaryDir>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;
using primerpicker::gui::test::table_rows;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

// From Molecule at 400 of TEST600, complement strand, saved with S
void complement_primer_a(MainWindow& window) {
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    child<QAction>(window, u"from_molecule").trigger();
    child<QCheckBox>(window, u"complement_box").setChecked(true);
    auto& field = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&field, u"400"_s);
    QTest::keyClick(&field, Qt::Key_Return);
    QTest::keyClick(&child<QWidget>(window, u"picking_view"), Qt::Key_S);
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

// Reference: PRIMER.exe, Match Active Primer for the primer at 400(C) of
// TEST600, target 100, the form's defaults (tests/legacy/complement/
// match_c400.PRN, recorded 2026-10-22): the form offers primers B on the
// top strand for a complement-strand primer A; one primer B, 73
TEST_CASE("Match Active Primer with a complement-strand primer A", "[gui][match][legacy]") {
    auto window = MainWindow{};
    complement_primer_a(window);
    child<QAction>(window, u"match_active_primer").trigger();
    CHECK_FALSE(child<QCheckBox>(window, u"match_complement_box").isChecked());  // the form's default
    auto& target = child<QLineEdit>(window, u"match_target_field");
    QTest::keyClicks(&target, u"100"_s);
    QTest::keyClick(&target, Qt::Key_Return);
    CHECK_FALSE(child<QDialog>(window, u"match_dialog").isVisible());
    CHECK(table_rows(window, u"match_table") ==
          QList<QList<QString>>{{u"73"_s, u"27"_s, u"50"_s, u"69"_s, u"2"_s, u"4"_s, u"1"_s, u"3"_s, u"2"_s, u"0"_s}});
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"REPORT"_s);
    REQUIRE_FALSE(window.print_report_to(path).has_value());
    CHECK(body(path) == body(legacy_file(u"complement/match_c400.PRN")));
}

// Reference: PRIMER.exe, the same primer A, target 450, "Complement
// strand (Y/N)" answered Y (complement/match_c400_compB.PRN): 7 primers B
// from 413 to 457
TEST_CASE("Match Active Primer with both primers on the complement strand", "[gui][match][legacy]") {
    auto window = MainWindow{};
    complement_primer_a(window);
    REQUIRE_FALSE(window.match_active_primer({.target = 450,
                                              .start_offset = 50,
                                              .end_offset = 50,
                                              .strand = primerpicker::sequence::Strand::complement})
                      .has_value());
    auto const rows = table_rows(window, u"match_table");
    REQUIRE(rows.size() == 7);
    CHECK(rows.front().front() == u"413"_s);
    CHECK(rows.back().front() == u"457"_s);
}
