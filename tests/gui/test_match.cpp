#include "gui/main_window.hpp"
#include "gui_primer_helpers.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QCheckBox>
#include <QDialog>
#include <QFile>
#include <QIODevice>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QString>
#include <QStringView>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;
using primerpicker::gui::test::direct_entry;
using primerpicker::gui::test::table_rows;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

// Primer/Create -> From Molecule at POSITION, top or complement strand
void from_molecule(MainWindow& window, QString const& position, bool const complement) {
    child<QAction>(window, u"from_molecule").trigger();
    child<QCheckBox>(window, u"complement_box").setChecked(complement);
    auto& field = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&field, position);
    QTest::keyClick(&field, Qt::Key_Return);
}

// the list of tests/legacy/search/match_active.tsv: position, bps from
// target, GC %, Tm, self-dimers (3', any), A:B dimers (3', any), runs,
// 3' G/C runs
auto recorded_list() -> QList<QList<QString>> {
    auto table = QFile{legacy_file(u"search/match_active.tsv")};
    REQUIRE(table.open(QIODevice::ReadOnly | QIODevice::Text));
    constexpr auto first_value = 8;
    auto rows = QList<QList<QString>>{};
    while (not table.atEnd()) {
        auto const row = QString::fromUtf8(table.readLine()).chopped(1);
        if (not row.startsWith(u'#') and not row.isEmpty()) {
            rows.append(row.split(u'\t').mid(first_value));
        }
    }
    return rows;
}

}  // namespace

// Reference: PRIMER.exe, Primer/Create -> Match Active Primer for the
// primer at 51 of TEST600.GB (top strand), target 400, the form's
// defaults (-50, +50): the Find Matching Primer list
// (tests/legacy/search/match_active.tsv, recorded 2026-10-18)
TEST_CASE("Match Active Primer lists the matching primers B", "[gui][match][legacy]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    from_molecule(window, u"51"_s, false);
    REQUIRE_FALSE(window.match_active_primer({.target = 400, .start_offset = 50, .end_offset = 50}).has_value());
    CHECK(child<QWidget>(window, u"match_view").isVisibleTo(&window));
    auto const recorded = recorded_list();
    REQUIRE(recorded.size() == 11);
    CHECK(table_rows(window, u"match_table") == recorded);
}

TEST_CASE("Match Active Primer asks for the target, with the settings' range", "[gui][match]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    from_molecule(window, u"51"_s, false);
    child<QAction>(window, u"match_active_primer").trigger();
    CHECK(child<QLineEdit>(window, u"match_start_field").text() == u"50"_s);
    CHECK(child<QLineEdit>(window, u"match_end_field").text() == u"50"_s);
    auto& target = child<QLineEdit>(window, u"match_target_field");
    QTest::keyClicks(&target, u"400"_s);
    QTest::keyClick(&target, Qt::Key_Return);
    CHECK_FALSE(child<QDialog>(window, u"match_dialog").isVisible());
    CHECK(table_rows(window, u"match_table").size() == 11);
}

TEST_CASE("Match Active Primer needs a top-strand primer made from the molecule", "[gui][match]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    direct_entry(window, u"CTCAAAGGCTCAGACATCAG"_s);
    CHECK(window.match_active_primer({.target = 400, .start_offset = 50, .end_offset = 50}).has_value());
    // a complement-strand primer A (primers B on the top strand):
    // recorded, not reproduced yet (match_variations/)
    from_molecule(window, u"400"_s, true);
    CHECK(window.match_active_primer({.target = 51, .start_offset = 50, .end_offset = 50}).has_value());
}

// Reference: "No Matching Primers Found" (tests/legacy/README.md, "No
// Matching Primers Found"): a target with no primer B in its window
TEST_CASE("Match Active Primer without a match shows PRIMER.exe's message", "[gui][match]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    from_molecule(window, u"51"_s, false);
    REQUIRE_FALSE(window.match_active_primer({.target = 400, .start_offset = 0, .end_offset = 0}).has_value());
    CHECK(child<QLabel>(window, u"match_message").text() == u"No Matching Primers Found"_s);
}

// Reference: PRIMER.exe, Match Active Primer for the primer at 50 of
// search/TRUN.GB, target 600, "Complement strand (Y/N)" answered N:
// rows 552 to 609, and the report printed to disk
// (tests/legacy/match_variations/TRUN_50_600_top_strand.PRN, recorded
// 2026-10-20)
TEST_CASE("Match Active Primer with primers B on the top strand", "[gui][match][legacy]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"search/TRUN.GB")).has_value());
    from_molecule(window, u"50"_s, false);
    child<QAction>(window, u"match_active_primer").trigger();
    auto& complement = child<QCheckBox>(window, u"match_complement_box");
    CHECK(complement.isChecked());  // the form's default
    complement.setChecked(false);
    auto& target = child<QLineEdit>(window, u"match_target_field");
    QTest::keyClicks(&target, u"600"_s);
    QTest::keyClick(&target, Qt::Key_Return);
    auto const rows = table_rows(window, u"match_table");
    REQUIRE_FALSE(rows.isEmpty());
    CHECK(rows.front().front() == u"552"_s);
    CHECK(rows.back().front() == u"609"_s);
}
