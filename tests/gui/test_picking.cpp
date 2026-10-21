#include "gui/colour_scheme.hpp"
#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QColor>
#include <QDir>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QPalette>
#include <QPushButton>
#include <QString>
#include <QStringView>
#include <QStyle>
#include <QTableView>
#include <QTest>
#include <QWidget>
#include <QtEnvironmentVariables>

#include <catch2/catch_test_macros.hpp>

#include <stdexcept>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::ColourScheme;
using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

enum struct Complement : bool { no, yes };

// Primer/Create -> From Molecule, with a molecule already loaded
void from_molecule(MainWindow& window, QString const& position, Complement const complement = Complement::no) {
    child<QAction>(window, u"from_molecule").trigger();
    auto& field = child<QLineEdit>(window, u"position_field");
    child<QCheckBox>(window, u"complement_box").setChecked(complement == Complement::yes);
    QTest::keyClicks(&field, position);
    QTest::keyClick(&field, Qt::Key_Return);
}

auto picking(MainWindow& window) -> QWidget& {
    return child<QWidget>(window, u"picking_view");
}

auto values(MainWindow const& window) -> QString {
    return child<QLabel>(window, u"picking_values").text();
}

auto verdict(MainWindow const& window) -> QString {
    return child<QLabel>(window, u"picking_verdict").text();
}

void key(MainWindow& window, Qt::Key const pressed, Qt::KeyboardModifiers const modifiers = Qt::NoModifier) {
    QTest::keyClick(&picking(window), pressed, modifiers);
}

auto list_origins(MainWindow const& window) -> QList<QString> {
    auto const& model = *child<QTableView>(window, u"primer_list").model();
    constexpr auto origin_column = 3;
    auto origins = QList<QString>{};
    for (auto row = 0; row < model.rowCount(); ++row) {
        origins.append(model.data(model.index(row, origin_column)).toString());
    }
    return origins;
}

void load_test600(MainWindow& window) {
    if (window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value()) {
        throw std::runtime_error{"cannot load TEST600.GB"};
    }
}

}  // namespace

// Reference: PRIMER.exe, From Molecule at 51 of TEST600, built-in
// criteria (tests/legacy/picking/picking.tsv, run explore_51)
TEST_CASE("From Molecule opens the picking view with the primer's values", "[gui][picking]") {
    auto window = MainWindow{};
    load_test600(window);
    from_molecule(window, u"51"_s);
    CHECK(picking(window).isVisibleTo(&window));
    CHECK(child<QLabel>(window, u"picking_header").text() == u"Primer length: 20   Molecule: TEST600"_s);
    CHECK(values(window) == u"Pos: 51     %GC: 50    Tm°C: 66    Hairpin: OK  Dimers: OK  Runs: !"_s);
    CHECK(verdict(window) == u"Meets criteria selected: No"_s);
}

// Same reference: Right gives 52; S puts it into the Primer List
// (save_52: "TEST600 ! 52 20 TCAAAGGCTCAGACATCAGG")
TEST_CASE("the arrows move the primer, S saves it to the Primer List", "[gui][picking]") {
    auto window = MainWindow{};
    load_test600(window);
    from_molecule(window, u"51"_s);
    key(window, Qt::Key_Right);
    CHECK(values(window) == u"Pos: 52     %GC: 50    Tm°C: 69    Hairpin: OK  Dimers: OK  Runs: !"_s);
    key(window, Qt::Key_S);
    CHECK(list_origins(window) == QList<QString>{u"TEST600 at 52"_s});
    CHECK(child<QLabel>(window, u"sequence_label").text() == u"5'-TCAAAGGCTCAGACATCAGG-3'"_s);
    CHECK(child<QWidget>(window, u"analysis_view").isVisibleTo(&window));
    CHECK_FALSE(picking(window).isVisibleTo(&window));
}

// Reference: PRIMER.exe saves with Esc too (esc_52)
TEST_CASE("Esc saves the primer too, as in PRIMER.exe", "[gui][picking]") {
    auto window = MainWindow{};
    load_test600(window);
    from_molecule(window, u"51"_s);
    key(window, Qt::Key_Right);
    key(window, Qt::Key_Escape);
    CHECK(list_origins(window) == QList<QString>{u"TEST600 at 52"_s});
}

// Reference: PRIMER.exe, Ctrl+Right from 51 gives 71, which meets the
// criteria (find_51)
TEST_CASE("Ctrl+Right finds the next primer that meets the criteria", "[gui][picking]") {
    auto window = MainWindow{};
    load_test600(window);
    from_molecule(window, u"51"_s);
    key(window, Qt::Key_Right, Qt::ControlModifier);
    CHECK(values(window) == u"Pos: 71     %GC: 55    Tm°C: 75    Hairpin: OK  Dimers: OK  Runs: OK"_s);
    CHECK(verdict(window) == u"Meets criteria selected: Yes"_s);
    child<QPushButton>(window, u"find_previous_button").click();
    CHECK(values(window).startsWith(u"Pos: 36 "_s));
}

// Reference: PRIMER.exe, from 559 of TEST600 Ctrl+Right finds nothing
// (find_end); + at 581 is refused (end_579)
TEST_CASE("the picking view shows PRIMER.exe's messages", "[gui][picking]") {
    auto window = MainWindow{};
    load_test600(window);
    from_molecule(window, u"559"_s);
    key(window, Qt::Key_Right, Qt::ControlModifier);
    CHECK(child<QLabel>(window, u"picking_message").text() == u"No Primers found"_s);
    CHECK(values(window).startsWith(u"Pos: 559 "_s));
    from_molecule(window, u"581"_s);
    key(window, Qt::Key_Plus);
    CHECK(child<QLabel>(window, u"picking_message").text() == u"Primer Too Long - Extends outside molecule"_s);
    CHECK(child<QLabel>(window, u"picking_header").text().startsWith(u"Primer length: 20 "_s));
    key(window, Qt::Key_Minus);
    CHECK(child<QLabel>(window, u"picking_message").text().isEmpty());
    CHECK(child<QLabel>(window, u"picking_header").text().startsWith(u"Primer length: 19 "_s));
}

// Reference: PRIMER.exe, complement strand at 400 (comp_400): Right
// gives 401; the screen says "*** Complementary Strand ***"
TEST_CASE("picking on the complement strand", "[gui][picking]") {
    auto window = MainWindow{};
    load_test600(window);
    from_molecule(window, u"400"_s, Complement::yes);
    CHECK(child<QLabel>(window, u"picking_strand").text() == u"*** Complementary Strand ***"_s);
    CHECK(values(window) == u"Pos: 400    %GC: 45!   Tm°C: 66    Hairpin: !   Dimers: OK  Runs: OK"_s);
    key(window, Qt::Key_Right);
    CHECK(values(window).startsWith(u"Pos: 401 "_s));
    key(window, Qt::Key_Return);
    CHECK(list_origins(window) == QList<QString>{u"TEST600 at 401 (complement)"_s});
}

// Reference: PRIMER.exe, Modify -> Position on a saved primer
// (modify_position: saved at 51, moved to 52, then 71, saved again)
TEST_CASE("Modify -> Position opens the picking view for the active primer", "[gui][picking]") {
    auto window = MainWindow{};
    load_test600(window);
    from_molecule(window, u"51"_s);
    key(window, Qt::Key_S);
    REQUIRE_FALSE(window.modify_position().has_value());
    CHECK(values(window).startsWith(u"Pos: 51 "_s));
    key(window, Qt::Key_Right);
    key(window, Qt::Key_Right, Qt::ControlModifier);
    key(window, Qt::Key_S);
    CHECK(list_origins(window) == QList<QString>{u"TEST600 at 51"_s, u"TEST600 at 71"_s});
}

TEST_CASE("Modify -> Position needs a primer made from the loaded molecule", "[gui][picking]") {
    auto window = MainWindow{};
    load_test600(window);
    CHECK(window.modify_position().has_value());
    REQUIRE_FALSE(window.load_primer(legacy_file(u"primers/M13EXT.PD2")).has_value());
    CHECK(window.modify_position().has_value());
}

// Not in PRIMER.exe (decided 2026-10-21): while the second primer is
// picked, the A:B dimers with the primer saved last. Reference values:
// PRIMER.exe's Primer Pair screen for TEST600 51 and 400(C), 3' dimers
// A:B 1, any 3 (tests/legacy/pairs/pair_screens.tsv)
TEST_CASE("while picking a second primer, the A:B dimers with the first are shown", "[gui][picking]") {
    auto window = MainWindow{};
    load_test600(window);
    from_molecule(window, u"51"_s);
    CHECK_FALSE(child<QLabel>(window, u"picking_dimers").isVisibleTo(&window));
    key(window, Qt::Key_S);
    from_molecule(window, u"400"_s, Complement::yes);
    auto const& dimers = child<QLabel>(window, u"picking_dimers");
    CHECK(dimers.isVisibleTo(&window));
    CHECK(dimers.text() ==
          u"A:B dimers with primer 1 (CTCAAAGGCTCAGACATCAG): 3' end matches 1 OK, adjacent bases 3 OK"_s);
}

namespace {

// the offscreen platform has no dark palette: a dark one is set by hand
auto dark_palette() -> QPalette {
    constexpr auto window_grey = 0x2d;
    constexpr auto base_grey = 0x1e;
    auto palette = QPalette{QColor{window_grey, window_grey, window_grey}};
    palette.setColor(QPalette::Base, QColor{base_grey, base_grey, base_grey});
    palette.setColor(QPalette::WindowText, Qt::white);
    palette.setColor(QPalette::Text, Qt::white);
    palette.setColor(QPalette::ButtonText, Qt::white);
    return palette;
}

void save_picking(QDir const& directory, ColourScheme const scheme) {
    primerpicker::gui::apply_colour_scheme(scheme);
    QApplication::setPalette(scheme == ColourScheme::dark ? dark_palette() : QApplication::style()->standardPalette());
    constexpr auto width = 900;
    constexpr auto height = 640;
    for (auto const meets : {false, true}) {
        auto window = MainWindow{};
        load_test600(window);
        window.resize(width, height);
        window.show();
        from_molecule(window, u"51"_s);
        if (meets) {
            key(window, Qt::Key_Right, Qt::ControlModifier);
        }
        QApplication::processEvents();
        auto const file = directory.filePath(
            u"picking_%1_%2.png"_s.arg(meets ? u"meets"_s : u"fails"_s, primerpicker::gui::colour_scheme_name(scheme)));
        CHECK(window.grab().save(file));
    }
}

}  // namespace

// Hidden test (not run by CTest), as tests/gui/test_screenshots.cpp:
//   PRIMERPICKER_SCREENSHOT_DIR=DIR build/debug/tests/gui/gui_tests '[.picking_screenshots]'
TEST_CASE("screenshots of the picking view", "[.picking_screenshots]") {
    auto const directory = qEnvironmentVariable("PRIMERPICKER_SCREENSHOT_DIR");
    REQUIRE_FALSE(directory.isEmpty());
    QApplication::setStyle(u"Fusion"_s);
    save_picking(QDir{directory}, ColourScheme::light);
    save_picking(QDir{directory}, ColourScheme::dark);
    primerpicker::gui::apply_colour_scheme(ColourScheme::system);
    QApplication::setPalette(QApplication::style()->standardPalette());
}
