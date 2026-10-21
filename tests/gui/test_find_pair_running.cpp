#include "gui/main_window.hpp"

#include <QString>
#include <QStringView>
#include <QTest>
#include <QThread>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

}  // namespace

// Found on 2026-10-21 (a rare failure of the Find Pair tests in release
// builds): the search thread may end before its result is shown, which
// happens when the window's event loop runs again. Until then the search
// must still count as running, or a caller waiting for its end finds no
// list.
TEST_CASE("Find Pair is running until its result is shown", "[gui][find pair]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    REQUIRE_FALSE(
        window.start_find_pair({.target_a = 150, .start_offset_a = 100, .target_b = 400, .end_offset_b = 100})
            .has_value());
    // the search of TEST600 takes milliseconds; no event is processed
    constexpr auto search_time_ms = 500;
    QThread::msleep(search_time_ms);
    CHECK(window.find_pair_running());
    constexpr auto timeout_ms = 10000;
    CHECK(QTest::qWaitFor([&window] { return not window.find_pair_running(); }, timeout_ms));
}
