#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QCheckBox>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>
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

enum struct Complement : bool { no, yes };

// Primer/Create -> From Molecule, with a molecule already loaded
void from_molecule(MainWindow& window, QString const& position, Complement const complement) {
    child<QAction>(window, u"from_molecule").trigger();
    auto& field = child<QLineEdit>(window, u"position_field");
    child<QCheckBox>(window, u"complement_box").setChecked(complement == Complement::yes);
    QTest::keyClicks(&field, position);
    QTest::keyClick(&field, Qt::Key_Return);
}

auto analysed_sequence(MainWindow const& window) -> QString {
    return child<QLabel>(window, u"sequence_label").text();
}

auto analysed_tm(MainWindow const& window) -> QString {
    auto const& model = *child<QTableView>(window, u"criteria_table").model();
    constexpr auto tm_row = 1;
    constexpr auto actual_column = 2;
    return model.data(model.index(tm_row, actual_column)).toString();
}

}  // namespace

TEST_CASE("a molecule file is loaded and named in the status bar", "[gui][molecule]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    CHECK(child<QLabel>(window, u"molecule_label").text() == u"Molecule: TEST600, 600 bp, linear"_s);
}

TEST_CASE("a file that cannot be read is not loaded", "[gui][molecule]") {
    auto window = MainWindow{};
    CHECK(window.load_molecule(legacy_file(u"no such file.gb")).has_value());
    CHECK(child<QLabel>(window, u"molecule_label").text().isEmpty());
}

// Reference: PRIMER.exe, Primer/Create -> From Molecule on TEST600.GB,
// built-in defaults (20 bases), Primer Pair screen
// (tests/legacy/pairs/pairs.tsv, recorded 2026-10-17): position 51 gives
// CTCAAAGGCTCAGACATCAG, Tm 66 °C
TEST_CASE("From Molecule makes the primer at a position of the top strand", "[gui][molecule]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    from_molecule(window, u"51"_s, Complement::no);
    CHECK_FALSE(child<QDialog>(window, u"from_molecule_dialog").isVisible());
    CHECK(analysed_sequence(window) == u"5'-CTCAAAGGCTCAGACATCAG-3'"_s);
    CHECK(analysed_tm(window) == u"66"_s);
}

// Same reference: position 400, complement strand, gives
// TGTATCCAGGTAATGCAAGG, Tm 66 °C
TEST_CASE("From Molecule makes the primer at a position of the complement strand", "[gui][molecule]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    from_molecule(window, u"400"_s, Complement::yes);
    CHECK(analysed_sequence(window) == u"5'-TGTATCCAGGTAATGCAAGG-3'"_s);
    CHECK(analysed_tm(window) == u"66"_s);
}

// Reference: PRIMER.exe refuses a primer made at 591 (20 bases in 600)
// with "Invalid Position", on circular and linear molecules
// (tests/legacy/README.md, "Circular molecules", recorded 2026-10-20)
TEST_CASE("From Molecule refuses a primer beyond the molecule's end", "[gui][molecule]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"circular/TCIRC.GB")).has_value());
    from_molecule(window, u"591"_s, Complement::no);
    CHECK(child<QDialog>(window, u"from_molecule_dialog").isVisible());
    CHECK(child<QLabel>(window, u"position_message").text() == u"Invalid Position"_s);
    CHECK(child<QLabel>(window, u"no_primer_label").isVisibleTo(&window));
}

TEST_CASE("the position field takes digits only", "[gui][molecule]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    child<QAction>(window, u"from_molecule").trigger();
    auto& field = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&field, u"4x0-0"_s);
    CHECK(field.text() == u"400"_s);
}
