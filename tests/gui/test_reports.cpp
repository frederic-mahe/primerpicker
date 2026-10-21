#include "core/analysis/legacy_homology.hpp"
#include "gui/main_window.hpp"
#include "gui_primer_helpers.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QByteArray>
#include <QCheckBox>
#include <QFile>
#include <QIODevice>
#include <QKeySequence>
#include <QLineEdit>
#include <QRegularExpression>
#include <QString>
#include <QStringView>
#include <QTemporaryDir>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::analysis::HomologyKind;
using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;
using primerpicker::gui::test::save_direct_entry;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

// a report without its dates (one per page header)
auto body(QString const& path) -> QByteArray {
    auto file = QFile{path};
    REQUIRE(file.open(QIODevice::ReadOnly));
    auto text = file.readAll();
    static auto const date = QRegularExpression{u"Report Date: [0-9/]+"_s};
    return QString::fromLatin1(text).replace(date, u"Report Date: DATE"_s).toLatin1();
}

// prints the report of the current screen into a temporary file, and
// compares it with the recording RECORDED, dates aside
void check_report(MainWindow& window, QStringView const recorded) {
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"REPORT"_s);
    REQUIRE_FALSE(window.print_report_to(path).has_value());
    CHECK(body(path) == body(legacy_file(recorded)));
}

void from_molecule(MainWindow& window, QString const& position) {
    child<QAction>(window, u"from_molecule").trigger();
    child<QCheckBox>(window, u"complement_box").setChecked(false);
    auto& field = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&field, position);
    QTest::keyClick(&field, Qt::Key_Return);
}

}  // namespace

TEST_CASE("Analyze -> Print to Disk prints the report of the current screen (Ctrl+P)", "[gui][report]") {
    auto const window = MainWindow{};
    CHECK(child<QAction>(window, u"print_report").shortcut() == QKeySequence{Qt::CTRL | Qt::Key_P});
}

// Reference: PRIMER.exe's reports printed to disk (tests/legacy/reports/
// reports.tsv, recorded 2026-10-19 and 2026-10-20), dates aside
TEST_CASE("Print to Disk: the Primer Pair report", "[gui][report][legacy]") {
    auto window = MainWindow{};
    save_direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    save_direct_entry(window, u"AACAGCTATGACCATG"_s);
    REQUIRE_FALSE(window.analyze_pair({.first = 0, .second = 1}).has_value());
    check_report(window, u"reports/PAIR.PRN");
}

TEST_CASE("Print to Disk: the Multiple Primers report", "[gui][report][legacy]") {
    auto window = MainWindow{};
    for (auto const& primer : {u"CTCAAAGGCTCAGACATCAG"_s, u"GATCCATGCAATGCTAGCTA"_s, u"ACGTTGCAGGTCAAGCTTGC"_s,
                               u"TTGGCCAAGTCGATCGGATC"_s}) {
        save_direct_entry(window, primer);
    }
    REQUIRE_FALSE(window.analyze_multi({std::size_t{0}, std::size_t{1}, std::size_t{2}, std::size_t{3}}).has_value());
    check_report(window, u"reports/MULTI.PRN");
}

TEST_CASE("Print to Disk: the Find Pair report", "[gui][report][legacy]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    REQUIRE_FALSE(window.find_pair({.target_a = 150, .start_offset_a = 100, .target_b = 400, .end_offset_b = 100})
                      .has_value());
    check_report(window, u"reports/FINDPAIR.PRN");
}

TEST_CASE("Print to Disk: the Uniqueness and False Priming reports", "[gui][report][legacy]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"homology/TREP2.GB")).has_value());
    from_molecule(window, u"51"_s);
    REQUIRE_FALSE(window.homology(HomologyKind::uniqueness).has_value());
    check_report(window, u"reports/UNIQUE.PRN");
    REQUIRE_FALSE(window.homology(HomologyKind::false_priming).has_value());
    check_report(window, u"reports/FALSEPRI.PRN");
}

TEST_CASE("Print to Disk: the Match Active Primer report", "[gui][report][legacy]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"search/TRUN.GB")).has_value());
    from_molecule(window, u"50"_s);
    REQUIRE_FALSE(window.match_active_primer({.target = 1020, .start_offset = 50, .end_offset = 50}).has_value());
    check_report(window, u"reports/match/TRUN_50_1020.PRN");
}

// Reference: PRIMER.exe prints nothing when there is no site ("Primer is
// unique": a box only; tests/legacy/README.md, "The Uniqueness and False
// Priming reports")
TEST_CASE("Print to Disk: nothing to print without a site", "[gui][report]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    from_molecule(window, u"51"_s);
    REQUIRE_FALSE(window.homology(HomologyKind::uniqueness).has_value());
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"REPORT"_s);
    CHECK(window.print_report_to(path).has_value());
    CHECK_FALSE(QFile::exists(path));
}
