#include "gui/main_window.hpp"

#include "gui_test_helpers.hpp"

#include <QAction>
#include <QCheckBox>
#include <QLabel>
#include <QLineEdit>
#include <QString>
#include <QStringView>
#include <QTest>
#include <QWidget>

#include <catch2/catch_test_macros.hpp>

#include <stdexcept>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

// From Molecule at 51 of TEST600 (saved with S), then the 5' extension
// GAATTC and a length of 20 bases in all, as PRIMER.exe's form makes it
// (51, length 20, extension GAATTC); then Modify -> Position
void picking_with_extension(MainWindow& window) {
    if (window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value()) {
        throw std::runtime_error{"cannot load TEST600.GB"};
    }
    child<QAction>(window, u"from_molecule").trigger();
    child<QCheckBox>(window, u"complement_box").setChecked(false);
    auto& position = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&position, u"51"_s);
    QTest::keyClick(&position, Qt::Key_Return);
    QTest::keyClick(&child<QWidget>(window, u"picking_view"), Qt::Key_S);
    child<QAction>(window, u"modify_extension").trigger();
    auto& extension = child<QLineEdit>(window, u"extension_field");
    QTest::keyClicks(&extension, u"GAATTC"_s);
    QTest::keyClick(&extension, Qt::Key_Return);
    REQUIRE_FALSE(window.change_length(20).has_value());
    REQUIRE_FALSE(window.modify_position().has_value());
}

auto values(MainWindow const& window) -> QString {
    return child<QLabel>(window, u"picking_values"_s).text();
}

}  // namespace

// Reference: PRIMER.exe, the picking screen of the primer made at 51 of
// TEST600 with the 5' extension GAATTC, 20 bases in all
// (tests/legacy/extension/picking_extension.tsv, recorded 2026-10-22)
TEST_CASE("the picking screen shows the values of a primer with its 5' extension", "[gui][picking]") {
    auto window = MainWindow{};
    picking_with_extension(window);
    CHECK(values(window) == u"Pos: 51     %GC: 45!   Tm°C: 47!   Hairpin: OK  Dimers: OK  Runs: !"_s);
    CHECK(child<QLabel>(window, u"picking_header"_s).text().contains(u"5'extension: GAATTC"_s));
    CHECK(child<QLabel>(window, u"picking_header"_s).text().contains(u"Primer length: 20"_s));
    QTest::keyClick(&child<QWidget>(window, u"picking_view"), Qt::Key_Right);
    QTest::keyClick(&child<QWidget>(window, u"picking_view"), Qt::Key_Plus);
    CHECK(values(window) == u"Pos: 52     %GC: 42!   Tm°C: 52!   Hairpin: OK  Dimers: OK  Runs: !"_s);
    QTest::keyClick(&child<QWidget>(window, u"picking_view"), Qt::Key_Minus);
    QTest::keyClick(&child<QWidget>(window, u"picking_view"), Qt::Key_Right, Qt::ControlModifier);
    CHECK(values(window) == u"Pos: 75     %GC: 50    Tm°C: 56    Hairpin: OK  Dimers: OK  Runs: OK"_s);
}
