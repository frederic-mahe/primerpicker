#include "gui/main_window.hpp"
#include "gui_primer_helpers.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QByteArray>
#include <QCheckBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFile>
#include <QIODevice>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QPushButton>
#include <QSpinBox>
#include <QString>
#include <QStringView>
#include <QTableView>
#include <QTemporaryDir>
#include <QTest>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cstddef>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;
using primerpicker::gui::test::direct_entry;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

auto read_all(QString const& path) -> QByteArray {
    auto file = QFile{path};
    REQUIRE(file.open(QIODevice::ReadOnly));
    return file.readAll();
}

auto gc_setting(MainWindow const& window) -> QString {
    auto const& model = *child<QTableView>(window, u"criteria_table").model();
    constexpr auto setting_column = 1;
    return model.data(model.index(0, setting_column)).toString();
}

void press_ok(MainWindow const& window, QStringView const buttons) {
    child<QDialogButtonBox>(window, buttons).button(QDialogButtonBox::Ok)->click();
}

}  // namespace

TEST_CASE("Modify -> Criteria changes the active primer's criteria only", "[gui][criteria]") {
    auto window = MainWindow{};
    direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    child<QAction>(window, u"modify_criteria").trigger();
    REQUIRE(child<QDialog>(window, u"criteria_dialog").isVisible());
    auto& gc_min = child<QSpinBox>(window, u"criteria_dialog_gc_min_box");
    CHECK(gc_min.value() == 50);  // the built-in default
    constexpr auto new_gc_min = 40;
    gc_min.setValue(new_gc_min);
    press_ok(window, u"criteria_dialog_buttons");
    CHECK_FALSE(child<QDialog>(window, u"criteria_dialog").isVisible());
    CHECK(gc_setting(window) == u"Min 40, Max 60"_s);
    // a new primer gets the System Defaults' criteria
    direct_entry(window, u"AACAGCTATGACCATG"_s);
    CHECK(gc_setting(window) == u"Min 50, Max 60"_s);
}

// Reference: PRIMER.exe's message for a value it refuses (from the
// strings of examples/PRIMER.exe)
TEST_CASE("Modify -> Criteria refuses a minimum above its maximum", "[gui][criteria]") {
    auto window = MainWindow{};
    direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    child<QAction>(window, u"modify_criteria").trigger();
    constexpr auto above_gc_max = 70;
    child<QSpinBox>(window, u"criteria_dialog_gc_min_box").setValue(above_gc_max);
    press_ok(window, u"criteria_dialog_buttons");
    CHECK(child<QDialog>(window, u"criteria_dialog").isVisible());
    CHECK(child<QLabel>(window, u"criteria_dialog_message").text() == u"Illegal value for defaults"_s);
    CHECK(gc_setting(window) == u"Min 50, Max 60"_s);
}

TEST_CASE("Modify -> Criteria without a primer asks for one", "[gui][criteria]") {
    auto window = MainWindow{};
    child<QAction>(window, u"modify_criteria").trigger();
    CHECK_FALSE(child<QDialog>(window, u"criteria_dialog").isVisible());
    CHECK(child<QLabel>(window, u"no_primer_label").isVisibleTo(&window));
}

// Reference: config/PD2_LENGTH25.CFG, saved by PRIMER.exe after the
// primer length was set to 25 in System Defaults (tests/legacy/README.md,
// "Settings file pd2.cfg")
TEST_CASE("Modify -> System Defaults: the primer length, saved as PRIMER.exe saves it", "[gui][settings][legacy]") {
    auto window = MainWindow{};
    child<QAction>(window, u"system_defaults").trigger();
    REQUIRE(child<QDialog>(window, u"system_defaults_dialog").isVisible());
    auto& length = child<QSpinBox>(window, u"primer_length_box");
    CHECK(length.value() == 20);
    constexpr auto new_length = 25;
    length.setValue(new_length);
    press_ok(window, u"system_defaults_buttons");
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    child<QAction>(window, u"from_molecule").trigger();
    auto& position = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&position, u"51"_s);
    QTest::keyClick(&position, Qt::Key_Return);
    CHECK(child<QLabel>(window, u"length_label").text() == u"25"_s);
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"PD2.CFG"_s);
    REQUIRE_FALSE(window.save_settings_to(path).has_value());
    CHECK(read_all(path) == read_all(legacy_file(u"config/PD2_LENGTH25.CFG")));
}

// Reference: config/PD2_WEIGHTS12345.CFG, the weights of the form set to
// 1 2 3 4 5 (GC mismatch, Tm mismatch, 3' dimers, any dimers, false
// priming; tests/legacy/README.md, "Settings file pd2.cfg")
TEST_CASE("Modify -> System Defaults: the Find Pair weights", "[gui][settings][legacy]") {
    auto window = MainWindow{};
    child<QAction>(window, u"system_defaults").trigger();
    constexpr auto weights = std::array{1, 2, 3, 4, 5};
    constexpr auto boxes = std::array<QStringView, 5>{u"weight_gc_box", u"weight_tm_box", u"weight_3prime_box",
                                                       u"weight_any_box", u"weight_false_priming_box"};
    for (auto index = std::size_t{0}; index < boxes.size(); ++index) {
        child<QSpinBox>(window, boxes.at(index)).setValue(weights.at(index));
    }
    press_ok(window, u"system_defaults_buttons");
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"PD2.CFG"_s);
    REQUIRE_FALSE(window.save_settings_to(path).has_value());
    CHECK(read_all(path) == read_all(legacy_file(u"config/PD2_WEIGHTS12345.CFG")));
}

// Reference: config/PD2_REPORT_NARROW.CFG (GC 40-65 %, Tm 50-70 °C,
// hairpin cutoff -3.5 kcal, dimers 4 and 6, runs 4, 3' G/C runs 2)
TEST_CASE("Modify -> System Defaults shows the loaded settings, and saves them unchanged", "[gui][settings][legacy]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_settings(legacy_file(u"config/PD2_REPORT_NARROW.CFG")).has_value());
    child<QAction>(window, u"system_defaults").trigger();
    CHECK(child<QSpinBox>(window, u"system_defaults_dialog_gc_min_box").value() == 40);
    CHECK(child<QSpinBox>(window, u"system_defaults_dialog_tm_max_box").value() == 70);
    constexpr auto narrow_cutoff = -3.5;
    constexpr auto tolerance = 1e-9;
    CHECK_THAT(child<QDoubleSpinBox>(window, u"system_defaults_dialog_hairpin_box").value(),
               Catch::Matchers::WithinAbs(narrow_cutoff, tolerance));
    CHECK(child<QSpinBox>(window, u"system_defaults_dialog_gc_run_box").value() == 2);
    press_ok(window, u"system_defaults_buttons");
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"PD2.CFG"_s);
    REQUIRE_FALSE(window.save_settings_to(path).has_value());
    CHECK(read_all(path) == read_all(legacy_file(u"config/PD2_REPORT_NARROW.CFG")));
}
