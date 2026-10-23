#include "gui/main_window.hpp"
#include "gui_primer_helpers.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QLabel>
#include <QString>
#include <QTableView>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;
using primerpicker::gui::test::save_direct_entry;

}  // namespace

// On macOS, item views do not emit activated() for Return (Return
// starts editing there); found by the Executables workflow on
// 2026-10-23. Blocking the view's signals does the same on every
// platform: Return and Enter must still pick the primer.
TEST_CASE("Return and Enter pick a primer of the list without the view's activated signal", "[gui][primer list]") {
    for (auto const key : {Qt::Key_Return, Qt::Key_Enter}) {
        CAPTURE(key);
        auto window = MainWindow{};
        save_direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
        save_direct_entry(window, u"AACAGCTATGACCATG"_s);
        auto& list = child<QTableView>(window, u"primer_list");
        list.setCurrentIndex(list.model()->index(0, 0));
        list.blockSignals(true);
        QTest::keyClick(&list, key);
        list.blockSignals(false);
        CHECK(child<QLabel>(window, u"sequence_label").text() == u"5'-GTAAAACGACGGCCAGT-3'"_s);
    }
}

TEST_CASE("the other keys still move in the Primer List", "[gui][primer list]") {
    auto window = MainWindow{};
    save_direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    save_direct_entry(window, u"AACAGCTATGACCATG"_s);
    auto& list = child<QTableView>(window, u"primer_list");
    list.setCurrentIndex(list.model()->index(0, 0));
    QTest::keyClick(&list, Qt::Key_Down);
    QTest::keyClick(&list, Qt::Key_Return);
    CHECK(child<QLabel>(window, u"sequence_label").text() == u"5'-AACAGCTATGACCATG-3'"_s);
}
