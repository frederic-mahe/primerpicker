#include "gui/main_window.hpp"
#include "gui/multi_primer_dialog.hpp"
#include "gui_primer_helpers.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFile>
#include <QIODevice>
#include <QList>
#include <QMap>
#include <QPushButton>
#include <QString>
#include <QStringView>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <optional>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::MultiSelection;
using primerpicker::gui::test::child;
using primerpicker::gui::test::save_direct_entry;
using primerpicker::gui::test::table_rows;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

auto read_rows(QStringView const file) -> QList<QList<QString>> {
    auto table = QFile{legacy_file(file)};
    REQUIRE(table.open(QIODevice::ReadOnly | QIODevice::Text));
    auto rows = QList<QList<QString>>{};
    while (not table.atEnd()) {
        auto const row = QString::fromUtf8(table.readLine()).chopped(1);
        if (not row.startsWith(u'#') and not row.isEmpty()) {
            rows.append(row.split(u'\t'));
        }
    }
    return rows;
}

// "1234", "1_2_": the list numbers typed in fields A to D ("_": empty)
auto selection_of(QString const& typed) -> MultiSelection {
    auto selection = MultiSelection{};
    for (auto field = std::size_t{0}; field < selection.size(); ++field) {
        auto const digit = typed.at(static_cast<qsizetype>(field));
        if (digit != u'_') {
            selection.at(field) = static_cast<std::size_t>(digit.digitValue() - 1);
        }
    }
    return selection;
}

auto field_index(QString const& field) -> int {
    return field.at(0).unicode() - u'A';
}

// A recorded screen: KEY "SELECTION LIST"; BALANCE: its rows (field,
// position, primer, length, % GC, Tm, 3' G/C runs, runs, hairpins);
// DIMERS: its cells (row, column, 3' end, any)
struct RecordedScreen {
    QString key;
    QList<QList<QString>> balance;
    QList<QList<QString>> dimers;
};

void check_balance(MainWindow const& window, QList<QList<QString>> const& balance) {
    auto expected = QList<QList<QString>>{};
    for (auto const& row : balance) {
        // the field, then from the primer on (the position is not shown)
        constexpr auto first_value = 2;
        expected.append(QList<QString>{row.at(0)} + row.mid(first_value));
    }
    CHECK(table_rows(window, u"multi_balance_table") == expected);
}

void check_dimers(MainWindow const& window, QList<QList<QString>> const& dimers) {
    auto const three_prime = table_rows(window, u"multi_3prime_table");
    auto const any = table_rows(window, u"multi_any_table");
    for (auto const& cell : dimers) {
        INFO(cell.at(0).toStdString() << cell.at(1).toStdString());
        auto const row = field_index(cell.at(0));
        auto const column = field_index(cell.at(1)) + 1;  // column 0 holds the row's field name
        CHECK(three_prime.at(row).at(column) == cell.at(2));
        CHECK(any.at(row).at(column) == cell.at(3));
        if (row != column - 1) {
            CHECK(three_prime.at(column - 1).at(row + 1) == u"."_s);  // below the diagonal
        }
    }
}

void check_screen(RecordedScreen const& screen) {
    INFO(screen.key.toStdString());
    auto const selection = screen.key.section(u' ', 0, 0);
    auto const primers = screen.key.section(u' ', 1, 1).split(u',');
    auto window = MainWindow{};
    for (auto const& primer : primers) {
        save_direct_entry(window, primer);
    }
    REQUIRE_FALSE(window.analyze_multi(selection_of(selection)).has_value());
    check_balance(window, screen.balance);
    check_dimers(window, screen.dimers);
}

}  // namespace

// Reference: PRIMER.exe, Analyze -> Multi Primer, built-in criteria:
// every recorded screen (tests/legacy/multi/multi_primer_balance.tsv and
// multi_primer_dimers.tsv, recorded 2026-10-19)
TEST_CASE("Multi Primer shows every recorded screen", "[gui][multi][legacy]") {
    constexpr auto key_fields = 2;
    auto balances = QMap<QString, QList<QList<QString>>>{};
    for (auto const& row : read_rows(u"multi/multi_primer_balance.tsv")) {
        balances[row.mid(0, key_fields).join(u' ')].append(row.mid(key_fields));
    }
    auto dimers = QMap<QString, QList<QList<QString>>>{};
    for (auto const& row : read_rows(u"multi/multi_primer_dimers.tsv")) {
        dimers[row.mid(0, key_fields).join(u' ')].append(row.mid(key_fields));
    }
    REQUIRE(balances.keys() == dimers.keys());
    for (auto const& key : balances.keys()) {
        check_screen({.key = key, .balance = balances.value(key), .dimers = dimers.value(key)});
    }
    CHECK(balances.size() == 9);
}

TEST_CASE("Analyze -> Multi Primer picks up to four primers from the Primer List", "[gui][multi]") {
    auto window = MainWindow{};
    save_direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    save_direct_entry(window, u"AACAGCTATGACCATG"_s);
    child<QAction>(window, u"analyze_multi_primer").trigger();
    auto& dialog = child<QDialog>(window, u"multi_primer_dialog");
    REQUIRE(dialog.isVisible());
    child<QComboBox>(window, u"multi_a_box").setCurrentIndex(1);  // index 0: empty
    child<QComboBox>(window, u"multi_b_box").setCurrentIndex(2);
    child<QDialogButtonBox>(window, u"multi_primer_buttons").button(QDialogButtonBox::Ok)->click();
    CHECK_FALSE(dialog.isVisible());
    CHECK(table_rows(window, u"multi_balance_table").size() == 2);
}

TEST_CASE("Multi Primer needs two filled fields", "[gui][multi]") {
    auto window = MainWindow{};
    save_direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    CHECK(window.analyze_multi({std::size_t{0}, std::nullopt, std::nullopt, std::nullopt}).has_value());
}
