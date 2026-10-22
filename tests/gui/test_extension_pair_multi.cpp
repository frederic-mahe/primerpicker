#include "gui/main_window.hpp"

#include "gui_primer_helpers.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QLineEdit>
#include <QList>
#include <QString>
#include <QStringView>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <optional>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;
using primerpicker::gui::test::direct_entry;
using primerpicker::gui::test::table_rows;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

struct ExtendedPrimer {
    QString primer;
    QString extension;
};

// Direct Entry of the primer, Modify -> 5'-extension, Save to List
void save_with_extension(MainWindow& window, ExtendedPrimer const& typed) {
    direct_entry(window, typed.primer);
    child<QAction>(window, u"modify_extension").trigger();
    auto& field = child<QLineEdit>(window, u"extension_field");
    QTest::keyClicks(&field, typed.extension);
    QTest::keyClick(&field, Qt::Key_Return);
    child<QAction>(window, u"save_to_list").trigger();
}

}  // namespace

// Reference: PRIMER.exe, Analyze -> Primer Pair for GAATTC +
// CTCAAAGGCTCAGA (made at 51 of TEST600, 20 bases in all) and
// TGTATCCAGGTAATGCAAGG (tests/legacy/extension/ext_pair.screens.txt,
// recorded 2026-10-22): the values of the whole primer, the Tm without
// the extension
TEST_CASE("Primer Pair with a 5' extension", "[gui][pair][legacy]") {
    auto window = MainWindow{};
    save_with_extension(window, {.primer = u"CTCAAAGGCTCAGA"_s, .extension = u"GAATTC"_s});
    primerpicker::gui::test::save_direct_entry(window, u"TGTATCCAGGTAATGCAAGG"_s);
    REQUIRE_FALSE(window.analyze_pair({.first = 0, .second = 1}).has_value());
    CHECK(table_rows(window, u"pair_table") ==
          QList<QList<QString>>{{u"Length"_s, u"20"_s, u"20"_s, u""_s},
                                {u"% GC"_s, u"45!"_s, u"45!"_s, u""_s},
                                {u"Tm °C"_s, u"47!"_s, u"66"_s, u"19°C difference"_s},
                                {u"Runs of G/C at 3' end"_s, u"0"_s, u"2"_s, u""_s},
                                {u"Runs of bases (any)"_s, u"3"_s, u"2"_s, u""_s},
                                {u"Hairpins"_s, u"OK"_s, u"!"_s, u"Hairpin(s) found"_s},
                                {u"3' end matches"_s, u"2"_s, u"2"_s, u"3!"_s},
                                {u"Any adj homol bases"_s, u"6"_s, u"4"_s, u"3"_s}});
}

// Reference: PRIMER.exe, Analyze -> Multi Primer for M13EXT.PD2 (GAATTC +
// GTAAAACGACGGCCAGT) and AACAGCTATGACCATG
// (tests/legacy/extension/ext_multi.screens.txt, recorded 2026-10-22)
TEST_CASE("Multi Primer with a 5' extension", "[gui][multi][legacy]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_primer(legacy_file(u"extension/P1.PD2")).has_value());
    child<QAction>(window, u"save_to_list").trigger();
    REQUIRE_FALSE(window.load_primer(legacy_file(u"extension/P2.PD2")).has_value());
    child<QAction>(window, u"save_to_list").trigger();
    REQUIRE_FALSE(window.analyze_multi({std::size_t{0}, std::size_t{1}, std::nullopt, std::nullopt}).has_value());
    auto const balance = table_rows(window, u"multi_balance_table");
    REQUIRE(balance.size() >= 2);
    CHECK(balance.at(0) == QList<QString>{u"A"_s, u"GAATTCGTAAAACGACGGCCAGT"_s, u"23"_s, u"47!"_s, u"64"_s, u"0"_s,
                                          u"4!"_s, u"OK"_s});
    CHECK(balance.at(1) ==
          QList<QString>{u"B"_s, u"AACAGCTATGACCATG"_s, u"16"_s, u"43!"_s, u"52!"_s, u"1"_s, u"2"_s, u"OK"_s});
    auto const three_prime = table_rows(window, u"multi_3prime_table");
    auto const any = table_rows(window, u"multi_any_table");
    CHECK(three_prime.at(0).at(1) == u"2"_s);  // A:A
    CHECK(any.at(0).at(1) == u"6"_s);
    CHECK(three_prime.at(0).at(2) == u"2"_s);  // A:B
    CHECK(any.at(0).at(2) == u"2"_s);
}
