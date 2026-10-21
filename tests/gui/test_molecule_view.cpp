#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QLabel>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QString>
#include <QStringView>
#include <QTest>
#include <QWidget>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

auto text_of(MainWindow const& window, QStringView const name) -> QString {
    return child<QPlainTextEdit>(window, name).toPlainText();
}

}  // namespace

// Reference: PRIMER.exe, Molecule-Ops -> View/Edit on pairs/TEST600.GB
// (tests/legacy/view_edit/, recorded 2026-10-21): the Sequence,
// Features and Annotations views
TEST_CASE("Molecule-Ops -> View/Edit shows the loaded molecule", "[gui][molecule]") {
    auto window = MainWindow{};
    CHECK(window.view_molecule() == u"Please load a molecule"_s);
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    REQUIRE_FALSE(window.view_molecule().has_value());
    CHECK(child<QWidget>(window, u"molecule_view").isVisibleTo(&window));
    CHECK(child<QLabel>(window, u"molecule_name_label").text() == u"TEST600"_s);
    CHECK(text_of(window, u"molecule_sequence").startsWith(
        u"      1  GGTCATTTTA AATGGACAAC CATCTATTGA GAATACTCAT GTCAACCTGG CTCAAAGGCT\n"
        u"     61  CAGACATCAG"_s));
    CHECK(text_of(window, u"molecule_features") == u"-- no features --"_s);
    CHECK(text_of(window, u"molecule_annotations").startsWith(u"LOCUS       TEST600      600 bp"_s));
}

// Same reference: Search for AAAT, "Number of Sites =  10", six sites a
// line (view_edit/find_aaat.screens.txt)
TEST_CASE("View/Edit: Find Sequence lists the sites on both strands", "[gui][molecule]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    REQUIRE_FALSE(window.view_molecule().has_value());
    auto& field = child<QLineEdit>(window, u"find_sequence_field");
    QTest::keyClicks(&field, u"aaat"_s);
    QTest::keyClick(&field, Qt::Key_Return);
    CHECK(child<QLabel>(window, u"find_sequence_count").text() == u"Number of Sites = 10"_s);
    CHECK(text_of(window, u"molecule_sites") ==
          u"        5C,        10 ,       113C,       190 ,       207C,       230C\n"
          u"      234C,       555 ,       582C,       589C"_s);
}

// Reference: PRIMER.exe, Molecule-Ops -> Save/Remove: "Molecule List is
// Empty" without a molecule; R unloads the molecule, the file stays
// (tests/legacy/view_edit/saveremove_*.screens.txt, recorded 2026-10-21)
TEST_CASE("Molecule-Ops -> Save/Remove: Remove unloads the molecule", "[gui][molecule]") {
    auto window = MainWindow{};
    CHECK(window.remove_molecule() == u"Molecule List is Empty"_s);
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    REQUIRE_FALSE(window.view_molecule().has_value());
    REQUIRE_FALSE(window.remove_molecule().has_value());
    CHECK(child<QLabel>(window, u"molecule_label").text().isEmpty());
    CHECK_FALSE(child<QWidget>(window, u"molecule_view").isVisibleTo(&window));
    CHECK(window.view_molecule() == u"Please load a molecule"_s);
}
