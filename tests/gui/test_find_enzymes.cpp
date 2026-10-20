#include "core/enzymes/site_search.hpp"
#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QCheckBox>
#include <QDialog>
#include <QFile>
#include <QIODevice>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QMap>
#include <QString>
#include <QStringView>
#include <QTableView>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::enzymes::Mismatches;
using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

auto table_rows(MainWindow const& window) -> QList<QList<QString>> {
    auto const& model = *child<QTableView>(window, u"sites_table").model();
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

// the recorded searches: "MOLECULE KIND ENZYME" -> rows (position,
// strand, match; no match for true sites)
auto recorded_searches() -> QMap<QString, QList<QList<QString>>> {
    auto table = QFile{legacy_file(u"enzymes/find_enzymes.tsv")};
    REQUIRE(table.open(QIODevice::ReadOnly | QIODevice::Text));
    auto searches = QMap<QString, QList<QList<QString>>>{};
    while (not table.atEnd()) {
        auto const row = QString::fromUtf8(table.readLine()).trimmed();
        if (row.startsWith(u'#') or row.isEmpty()) {
            continue;
        }
        auto const fields = row.split(u'\t');
        // molecule, kind, enzyme, recognition sequence, rank, position,
        // strand, match
        constexpr auto field_count = 8;
        constexpr auto position_field = 5;
        constexpr auto strand_field = 6;
        constexpr auto match_field = 7;
        REQUIRE(fields.size() == field_count);
        auto const key = fields.at(0) + u' ' + fields.at(1) + u' ' + fields.at(2);
        auto const match = fields.at(match_field) == u"-"_s ? QString{} : fields.at(match_field);
        searches[key].append({fields.at(position_field), fields.at(strand_field), match});
    }
    return searches;
}

// KEY: "MOLECULE KIND ENZYME"
void check_search(QString const& key, QList<QList<QString>> const& rows) {
    INFO(key.toStdString());
    auto const fields = key.split(u' ');
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"enzymes/" + fields.at(0))).has_value());
    auto const mismatches = fields.at(1) == u"almost"_s ? Mismatches::at_most_one : Mismatches::none;
    REQUIRE_FALSE(window.find_enzyme(fields.at(2), mismatches).has_value());
    CHECK(child<QWidget>(window, u"sites_view").isVisibleTo(&window));
    CHECK(table_rows(window) == rows);
    CHECK(child<QLabel>(window, u"sites_header").text().endsWith(u"Number of Sites = %1"_s.arg(rows.size())));
}

}  // namespace

// Reference: PRIMER.exe, Molecule-Ops -> Find Enzymes, True Sites and
// Almost Sites, on linear and circular molecules: every recorded list
// (tests/legacy/enzymes/find_enzymes.tsv, recorded 2026-10-18)
TEST_CASE("Find Enzymes lists the sites of every recorded search", "[gui][enzymes][legacy]") {
    auto const searches = recorded_searches();
    for (auto const& [key, rows] : searches.asKeyValueRange()) {
        check_search(key, rows);
    }
    CHECK(searches.size() == 12);
}

// Same reference, header of the Almost Sites screen: the reverse
// complement of a non-palindromic site in parentheses (as the CLI's
// find-enzymes, tests/cli/test_find_enzymes.sh)
TEST_CASE("the header of the Almost Sites list", "[gui][enzymes]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"enzymes/TENZ.GB")).has_value());
    REQUIRE_FALSE(window.find_enzyme(u"AlwI"_s, Mismatches::at_most_one).has_value());
    CHECK(child<QLabel>(window, u"sites_header").text() == u"Search for: AlwI GGATC (GATCC)    Number of Sites = 18"_s);
}

TEST_CASE("an unknown enzyme is refused", "[gui][enzymes]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"enzymes/TENZ.GB")).has_value());
    CHECK(window.find_enzyme(u"NoSuchI"_s, Mismatches::none).has_value());
}

TEST_CASE("Molecule-Ops -> Find Enzymes asks for the enzyme and the kind of sites", "[gui][enzymes]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"enzymes/TENZ.GB")).has_value());
    child<QAction>(window, u"find_enzymes").trigger();
    auto& field = child<QLineEdit>(window, u"enzyme_field");
    child<QCheckBox>(window, u"almost_box").setChecked(true);
    QTest::keyClicks(&field, u"ApaI"_s);
    QTest::keyClick(&field, Qt::Key_Return);
    CHECK_FALSE(child<QDialog>(window, u"find_enzymes_dialog").isVisible());
    CHECK(table_rows(window).size() == 3);
}

TEST_CASE("an unknown enzyme keeps the dialog open with a message", "[gui][enzymes]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"enzymes/TENZ.GB")).has_value());
    child<QAction>(window, u"find_enzymes").trigger();
    auto& field = child<QLineEdit>(window, u"enzyme_field");
    QTest::keyClicks(&field, u"NoSuchI"_s);
    QTest::keyClick(&field, Qt::Key_Return);
    CHECK(child<QDialog>(window, u"find_enzymes_dialog").isVisible());
    CHECK_FALSE(child<QLabel>(window, u"enzyme_message").text().isEmpty());
}
