#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QFile>
#include <QIODevice>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QString>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

// a row of tests/legacy/inspect/inspect.tsv: primer, the counts, then
// the hairpin energy and the three lines of its drawing ("-": none)
constexpr auto primer_field = 0;
constexpr auto hairpin_field = 6;
constexpr auto top_field = 7;
constexpr auto marks_field = 8;
constexpr auto bottom_field = 9;
constexpr auto field_count = 10;

void direct_entry(MainWindow& window, QString const& primer) {
    child<QAction>(window, u"direct_entry").trigger();
    auto& field = child<QLineEdit>(window, u"direct_entry_field");
    QTest::keyClicks(&field, primer);
    QTest::keyClick(&field, Qt::Key_Return);
}

}  // namespace

// Reference: PRIMER.exe, Direct Entry -> Inspect -> H (Hairpins
// screen), built-in defaults: every recorded primer of
// tests/legacy/inspect/inspect.tsv (recorded 2026-10-17); the drawing
// keeps the screen's columns
TEST_CASE("the hairpin drawing of every recorded primer", "[gui][hairpin][legacy]") {
    auto table = QFile{QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u"/inspect/inspect.tsv"_s};
    REQUIRE(table.open(QIODevice::ReadOnly | QIODevice::Text));
    auto window = MainWindow{};
    auto checked = std::size_t{0};
    while (not table.atEnd()) {
        auto const row = QString::fromUtf8(table.readLine()).chopped(1);
        if (row.startsWith(u'#')) {
            continue;
        }
        auto const fields = row.split(u'\t');
        REQUIRE(fields.size() == field_count);
        auto const expected = fields.at(hairpin_field) == u"-"_s
                                  ? u"No Hairpins Found"_s
                                  : QList<QString>{fields.at(top_field), fields.at(marks_field),
                                                   fields.at(bottom_field)}
                                        .join(u'\n');
        INFO(fields.at(primer_field).toStdString());
        direct_entry(window, fields.at(primer_field));
        CHECK(child<QLabel>(window, u"hairpin_drawing").text() == expected);
        ++checked;
    }
    CHECK(checked == 111);
}
