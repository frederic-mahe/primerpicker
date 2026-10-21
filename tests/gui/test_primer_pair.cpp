#include "gui/main_window.hpp"
#include "gui_primer_helpers.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFile>
#include <QIODevice>
#include <QList>
#include <QPushButton>
#include <QString>
#include <QStringView>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;
using primerpicker::gui::test::save_direct_entry;
using primerpicker::gui::test::table_rows;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

// "-" in the recordings: an empty comment
auto shown(QString const& text) -> QString {
    return text == u"-"_s ? QString{} : text;
}

// A row of pairs/pair_screens*.tsv: primer A, primer B, then three
// values (A, B, comment) for length, % GC, Tm, 3' G/C runs, runs and
// hairpins, then A:A, B:B and A:B for 3' end matches and for adjacent
// bases
void check_screen(QList<QString> const& fields) {
    constexpr auto field_count = 26;
    REQUIRE(fields.size() == field_count);
    INFO(fields.at(0).toStdString() << ' ' << fields.at(1).toStdString());
    auto window = MainWindow{};
    save_direct_entry(window, fields.at(0));
    save_direct_entry(window, fields.at(1));
    REQUIRE_FALSE(window.analyze_pair({.first = 0, .second = 1}).has_value());
    auto expected = QList<QList<QString>>{};
    auto const labels = QList<QString>{u"Length"_s,
                                       u"% GC"_s,
                                       u"Tm °C"_s,
                                       u"Runs of G/C at 3' end"_s,
                                       u"Runs of bases (any)"_s,
                                       u"Hairpins"_s,
                                       u"3' end matches"_s,
                                       u"Any adj homol bases"_s};
    constexpr auto first_value = 2;
    constexpr auto three_value_rows = 6;
    constexpr auto values_per_row = 3;  // A, B, comment; or A:A, B:B, A:B
    for (auto row = 0; row < labels.size(); ++row) {
        auto const start = first_value + (values_per_row * row);
        auto const third = row < three_value_rows ? shown(fields.at(start + 2)) : fields.at(start + 2);
        expected.append({labels.at(row), fields.at(start), fields.at(start + 1), third});
    }
    CHECK(table_rows(window, u"pair_table") == expected);
}

auto recorded_screens(QStringView const file) -> QList<QList<QString>> {
    auto table = QFile{legacy_file(file)};
    REQUIRE(table.open(QIODevice::ReadOnly | QIODevice::Text));
    auto screens = QList<QList<QString>>{};
    while (not table.atEnd()) {
        auto const row = QString::fromUtf8(table.readLine()).chopped(1);
        if (not row.startsWith(u'#') and not row.isEmpty()) {
            screens.append(row.split(u'\t'));
        }
    }
    return screens;
}

}  // namespace

// Reference: PRIMER.exe, Analyze -> Primer Pair of two primers typed with
// Direct Entry or made from molecules, built-in criteria: every recorded
// screen (tests/legacy/pairs/pair_screens.tsv, recorded 2026-10-18, and
// pair_screens_2.tsv)
TEST_CASE("Primer Pair shows every recorded screen", "[gui][pair][legacy]") {
    auto const screens = recorded_screens(u"pairs/pair_screens.tsv") + recorded_screens(u"pairs/pair_screens_2.tsv");
    for (auto const& screen : screens) {
        check_screen(screen);
    }
    CHECK(screens.size() == 42);
}

TEST_CASE("Analyze -> Primer Pair picks A and B from the Primer List", "[gui][pair]") {
    auto window = MainWindow{};
    save_direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    save_direct_entry(window, u"AACAGCTATGACCATG"_s);
    save_direct_entry(window, u"CTCAAAGGCTCAGACATCAG"_s);
    child<QAction>(window, u"analyze_primer_pair").trigger();
    auto& dialog = child<QDialog>(window, u"primer_pair_dialog");
    REQUIRE(dialog.isVisible());
    child<QComboBox>(window, u"pair_a_box").setCurrentIndex(2);
    child<QComboBox>(window, u"pair_b_box").setCurrentIndex(0);
    child<QDialogButtonBox>(window, u"primer_pair_buttons").button(QDialogButtonBox::Ok)->click();
    CHECK_FALSE(dialog.isVisible());
    CHECK(table_rows(window, u"pair_table").at(0) ==
          QList<QString>{u"Length"_s, u"20"_s, u"17"_s, u"3 base difference"_s});
}

TEST_CASE("Primer Pair needs two primers in the list", "[gui][pair]") {
    auto window = MainWindow{};
    save_direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    CHECK(window.analyze_pair({.first = 0, .second = 1}).has_value());
}
