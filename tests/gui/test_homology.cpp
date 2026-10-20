#include "core/analysis/legacy_homology.hpp"
#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QCheckBox>
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

using primerpicker::analysis::HomologyKind;
using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

auto table_rows(MainWindow const& window) -> QList<QList<QString>> {
    auto const& model = *child<QTableView>(window, u"homology_table").model();
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

// the recorded screens: "MOLECULE POSITION STRAND KIND" -> sites
// (homologous bases, position, strand, Tm); "-" rows: no site
auto recorded_screens() -> QMap<QString, QList<QList<QString>>> {
    auto table = QFile{legacy_file(u"homology/homology.tsv")};
    REQUIRE(table.open(QIODevice::ReadOnly | QIODevice::Text));
    // molecule, primer, position, strand, kind, homologous bases, site
    // position, site strand, Tm
    constexpr auto field_count = 9;
    constexpr auto position_field = 2;
    constexpr auto strand_field = 3;
    constexpr auto kind_field = 4;
    constexpr auto first_site_field = 5;
    auto screens = QMap<QString, QList<QList<QString>>>{};
    while (not table.atEnd()) {
        auto const row = QString::fromUtf8(table.readLine()).trimmed();
        if (row.startsWith(u'#') or row.isEmpty()) {
            continue;
        }
        auto const fields = row.split(u'\t');
        REQUIRE(fields.size() == field_count);
        auto const key = QList<QString>{fields.at(0), fields.at(position_field), fields.at(strand_field),
                                        fields.at(kind_field)}
                             .join(u' ');
        if (fields.at(first_site_field) == u"-"_s) {
            screens.insert(key, {});  // a screen without site
            continue;
        }
        screens[key].append(fields.mid(first_site_field));
    }
    return screens;
}

// KEY: "MOLECULE POSITION STRAND KIND"
void check_screen(QString const& key, QList<QList<QString>> const& sites) {
    INFO(key.toStdString());
    auto const fields = key.split(u' ');
    auto window = MainWindow{};
    auto const directory = fields.at(0).startsWith(u"TREP"_s) ? u"homology/"_s : u"pairs/"_s;
    REQUIRE_FALSE(window.load_molecule(legacy_file(directory + fields.at(0))).has_value());
    child<QAction>(window, u"from_molecule").trigger();
    child<QCheckBox>(window, u"complement_box").setChecked(fields.at(2) == u"complement"_s);
    auto& position = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&position, fields.at(1));
    QTest::keyClick(&position, Qt::Key_Return);
    auto const kind = fields.at(3) == u"uniqueness"_s ? HomologyKind::uniqueness : HomologyKind::false_priming;
    REQUIRE_FALSE(window.homology(kind).has_value());
    CHECK(child<QWidget>(window, u"homology_view").isVisibleTo(&window));
    CHECK(table_rows(window) == sites);
    auto const message = child<QLabel>(window, u"homology_message").text();
    CHECK(message.isEmpty() == not sites.isEmpty());
}

}  // namespace

// Reference: PRIMER.exe, Uniqueness (U) and False Priming (F) screens of
// primers made from molecules, built-in settings: every recorded screen
// (tests/legacy/homology/homology.tsv, recorded 2026-10-18)
TEST_CASE("Uniqueness and False Priming of every recorded primer", "[gui][homology][legacy]") {
    auto const screens = recorded_screens();
    for (auto const& [key, sites] : screens.asKeyValueRange()) {
        check_screen(key, sites);
    }
    CHECK(screens.size() == 7);
}

// Same reference: "Primer is unique" (TEST600.GB, 51)
TEST_CASE("a unique primer shows PRIMER.exe's message", "[gui][homology]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    child<QAction>(window, u"from_molecule").trigger();
    auto& position = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&position, u"51"_s);
    QTest::keyClick(&position, Qt::Key_Return);
    REQUIRE_FALSE(window.homology(HomologyKind::uniqueness).has_value());
    CHECK(child<QLabel>(window, u"homology_message").text() == u"Primer is unique"_s);
}

TEST_CASE("Uniqueness needs a primer made from the loaded molecule", "[gui][homology]") {
    auto window = MainWindow{};
    child<QAction>(window, u"direct_entry").trigger();
    auto& field = child<QLineEdit>(window, u"direct_entry_field");
    QTest::keyClicks(&field, u"CTCAAAGGCTCAGACATCAG"_s);
    QTest::keyClick(&field, Qt::Key_Return);
    CHECK(window.homology(HomologyKind::uniqueness).has_value());
}
