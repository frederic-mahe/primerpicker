#include "gui/main_window.hpp"

#include "gui_test_helpers.hpp"

#include <QAction>
#include <QCheckBox>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>
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

// the fields of PRIMER.exe's Create Primer form
struct Form {
    QString position;
    QString length;
    QString extension;
};

// Primer/Create -> From Molecule on TEST600, top strand
void from_molecule(MainWindow& window, Form const& form) {
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    child<QAction>(window, u"from_molecule").trigger();
    child<QCheckBox>(window, u"complement_box").setChecked(false);
    auto& length = child<QLineEdit>(window, u"length_field");
    length.clear();
    QTest::keyClicks(&length, form.length);
    QTest::keyClicks(&child<QLineEdit>(window, u"extension_field_form"), form.extension);
    auto& position = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&position, form.position);
    QTest::keyClick(&position, Qt::Key_Return);
}

}  // namespace

// Reference: PRIMER.exe, the Create Primer form (tests/legacy/extension/
// ext_picking3.screens.txt, recorded 2026-10-22): "Primer length (incl
// 5' ext): 20" (the System Defaults' length), "5' extension to be
// added:" empty
TEST_CASE("From Molecule offers the primer length and a 5' extension", "[gui][molecule]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    child<QAction>(window, u"from_molecule").trigger();
    CHECK(child<QLineEdit>(window, u"length_field").text() == u"20"_s);
    CHECK(child<QLineEdit>(window, u"extension_field_form").text().isEmpty());
}

// Reference: PRIMER.exe, From Molecule at 51 of TEST600, length 20,
// extension GAATTC (tests/legacy/extension/picking_extension.tsv,
// recorded 2026-10-22): the picking screen of GAATTC + 14 bases
TEST_CASE("From Molecule makes a primer with a 5' extension, its length included", "[gui][molecule][legacy]") {
    auto window = MainWindow{};
    from_molecule(window, {.position = u"51"_s, .length = u"20"_s, .extension = u"GAATTC"_s});
    CHECK_FALSE(child<QDialog>(window, u"from_molecule_dialog").isVisible());
    CHECK(child<QLabel>(window, u"picking_values").text() ==
          u"Pos: 51     %GC: 45!   Tm°C: 47!   Hairpin: OK  Dimers: OK  Runs: !"_s);
    CHECK(child<QLabel>(window, u"picking_header").text().contains(u"5'extension: GAATTC"_s));
    // Modify -> Length 25 then gives GAATTC + 19 bases (extension/ext_length2)
    QTest::keyClick(&child<QWidget>(window, u"picking_view"), Qt::Key_S);
    REQUIRE_FALSE(window.change_length(25).has_value());
    CHECK(child<QLabel>(window, u"sequence_label").text() == u"5'-GAATTC-CTCAAAGGCTCAGACATCA-3'"_s);
}

// Reference: PRIMER.exe, From Molecule at 51 of TEST600 with the form's
// length 25 (tests/legacy/modify_length/length_51_25: the primer of 25
// bases at 51 is CTCAAAGGCTCAGACATCAGGGCTC)
TEST_CASE("From Molecule makes a primer of the length typed", "[gui][molecule][legacy]") {
    auto window = MainWindow{};
    from_molecule(window, {.position = u"51"_s, .length = u"25"_s, .extension = {}});
    QTest::keyClick(&child<QWidget>(window, u"picking_view"), Qt::Key_S);
    CHECK(child<QLabel>(window, u"sequence_label").text() == u"5'-CTCAAAGGCTCAGACATCAGGGCTC-3'"_s);
}
