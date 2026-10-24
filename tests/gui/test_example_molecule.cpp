#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QLabel>
#include <QLineEdit>
#include <QString>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

}  // namespace

// The example molecule is pUC19 (data/molecules/pUC19.gb, GenBank
// L09137.2): 2686 bp, circular; the M13 forward (-20) primer
// GTAAAACGACGGCCAGT starts at base 379 of the record.
TEST_CASE("Help -> Example Molecule loads pUC19", "[gui][example]") {
    auto window = MainWindow{};
    child<QAction>(window, u"example_molecule").trigger();
    auto const label = child<QLabel>(window, u"molecule_label").text();
    CHECK(label.contains(u"2686 bp"_s));
    CHECK(label.contains(u"circular"_s));
}

TEST_CASE("primers can be made from the example molecule", "[gui][example]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_example_molecule().has_value());
    child<QAction>(window, u"from_molecule").trigger();
    auto& position = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&position, u"379"_s);
    QTest::keyClick(&position, Qt::Key_Return);
    CHECK(child<QLabel>(window, u"sequence_label").text().startsWith(u"5'-GTAAAACGACGGCCAGT"_s));
}
