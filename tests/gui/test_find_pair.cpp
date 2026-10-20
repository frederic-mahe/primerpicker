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
#include <QTest>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

auto table_rows(MainWindow const& window) -> QList<QList<QString>> {
    auto const& model = *child<QTableView>(window, u"pairs_table").model();
    auto rows = QList<QList<QString>>{};
    for (auto row = 0; row < model.rowCount(); ++row) {
        auto texts = QList<QString>{};
        for (auto column = 0; column < model.columnCount(); ++column) {
            texts.append(model.data(model.index(row, column)).toString());
        }
        rows.append(texts);
    }
    return rows;
}

// the ranked list of tests/legacy/search/find_pair.tsv: position A and
// B, GC A and B, Tm A and B, max 3' dimers, max dimers, false priming
// Tm A and B
auto recorded_list() -> QList<QList<QString>> {
    auto table = QFile{legacy_file(u"search/find_pair.tsv")};
    REQUIRE(table.open(QIODevice::ReadOnly | QIODevice::Text));
    constexpr auto first_value_field = 7;
    auto rows = QList<QList<QString>>{};
    while (not table.atEnd()) {
        auto const row = QString::fromUtf8(table.readLine()).trimmed();
        if (row.startsWith(u'#') or row.isEmpty()) {
            continue;
        }
        rows.append(row.split(u'\t').mid(first_value_field));
    }
    return rows;
}

auto label_text(MainWindow const& window, QStringView const name) -> QString {
    return child<QLabel>(window, name).text();
}

}  // namespace

// Reference: PRIMER.exe, Primer/Create -> Find Pair on TEST600.GB,
// targets 150 and 400, built-in defaults: the ranked list
// (tests/legacy/search/find_pair.tsv, recorded 2026-10-18)
TEST_CASE("Find Pair ranks the pairs as PRIMER.exe", "[gui][find pair][legacy]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    REQUIRE_FALSE(window.find_pair({.target_a = 150, .start_offset_a = 100, .target_b = 400, .end_offset_b = 100})
                      .has_value());
    CHECK(child<QWidget>(window, u"find_pair_view").isVisibleTo(&window));
    CHECK(table_rows(window) == recorded_list());
}

// Reference: "No Matching Primers Found" on TLONG.GB, targets 120 and
// 300 (tests/legacy/search/find_pair_no_match.tsv, recorded 2026-10-19)
TEST_CASE("Find Pair without a pair shows PRIMER.exe's message", "[gui][find pair][legacy]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"enzymes/TLONG.GB")).has_value());
    REQUIRE_FALSE(window.find_pair({.target_a = 120, .start_offset_a = 100, .target_b = 300, .end_offset_b = 100})
                      .has_value());
    CHECK(label_text(window, u"find_pair_message") == u"No Matching Primers Found"_s);
    CHECK(table_rows(window).isEmpty());
}

// Reference: PD2_RELAXED.CFG, TLONG.GB, targets 400 and 800: offsets
// 40 and 40 give 41 and 41 possible primers and 626 pairs; offsets 99
// and 5, "More than 100 sites found"
// (tests/legacy/search/find_pair_limits.tsv, recorded 2026-10-19)
TEST_CASE("Find Pair counts the possible primers and pairs", "[gui][find pair][legacy]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_settings(legacy_file(u"config/PD2_RELAXED.CFG")).has_value());
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"enzymes/TLONG.GB")).has_value());
    REQUIRE_FALSE(window.find_pair({.target_a = 400, .start_offset_a = 40, .target_b = 800, .end_offset_b = 40})
                      .has_value());
    CHECK(label_text(window, u"find_pair_counts") ==
          u"Primer A: 41 possible primers    Primer B: 41 possible primers    Pairs: 626"_s);
    CHECK(table_rows(window).size() == 626);
    REQUIRE_FALSE(window.find_pair({.target_a = 400, .start_offset_a = 99, .target_b = 800, .end_offset_b = 5})
                      .has_value());
    CHECK(label_text(window, u"find_pair_message") ==
          u"More than 100 sites found.  Suggest you redefine matching limits."_s);
    CHECK(table_rows(window).isEmpty());
}

TEST_CASE("Primer/Create -> Find Pair asks for the targets, with the settings' offsets", "[gui][find pair]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    child<QAction>(window, u"find_pair").trigger();
    CHECK(child<QLineEdit>(window, u"start_offset_field").text() == u"100"_s);
    CHECK(child<QLineEdit>(window, u"end_offset_field").text() == u"100"_s);
    QTest::keyClicks(&child<QLineEdit>(window, u"target_a_field"), u"150"_s);
    auto& target_b = child<QLineEdit>(window, u"target_b_field");
    QTest::keyClicks(&target_b, u"400"_s);
    QTest::keyClick(&target_b, Qt::Key_Return);
    CHECK_FALSE(child<QDialog>(window, u"find_pair_dialog").isVisible());
    CHECK(table_rows(window).size() == 10);
}
