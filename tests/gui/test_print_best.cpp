#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QApplication>
#include <QByteArray>
#include <QFile>
#include <QIODevice>
#include <QInputDialog>
#include <QRegularExpression>
#include <QString>
#include <QStringView>
#include <QTemporaryDir>
#include <QTimer>
#include <QWidget>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

// a report without its dates
auto body(QString const& path) -> QByteArray {
    auto file = QFile{path};
    REQUIRE(file.open(QIODevice::ReadOnly));
    auto text = file.readAll();
    static auto const date = QRegularExpression{u"Report Date: [0-9/]+"_s};
    return QString::fromLatin1(text).replace(date, u"Report Date: DATE"_s).toLatin1();
}

// Find Pair on search/TFP.GB, targets 150 and 400: 30 pairs
void find_tfp_pairs(MainWindow& window) {
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"search/TFP.GB")).has_value());
    REQUIRE_FALSE(window.find_pair({.target_a = 150, .start_offset_a = 100, .target_b = 400, .end_offset_b = 100})
                      .has_value());
}

}  // namespace

// Reference: PRIMER.exe, Find Pair on search/TFP.GB, targets 150 and
// 400, Print to Disk with "Print best: [25]" and "[9999]"
// (tests/legacy/reports/find_pair_best/, recorded 2026-10-21)
TEST_CASE("Print to Disk on the Find Pair list prints the best N pairs", "[gui][report][legacy]") {
    auto window = MainWindow{};
    find_tfp_pairs(window);
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"REPORT"_s);
    for (auto const best : {25, 9999}) {
        INFO(best);
        REQUIRE_FALSE(window.print_report_to(path, static_cast<std::size_t>(best)).has_value());
        CHECK(body(path) ==
              body(legacy_file(u"reports/find_pair_best/TFP_150_400_best%1.PRN"_s.arg(best))));
    }
}

// Reference: PRIMER.exe's box "30 Pairs found.  Print best: [20  ]"
// (tests/legacy/README.md, "Reports of more than 20 pairs")
TEST_CASE("Print to Disk on the Find Pair list asks how many pairs to print", "[gui][report]") {
    auto window = MainWindow{};
    find_tfp_pairs(window);
    auto label = QString{};
    auto value = 0;
    auto maximum = 0;
    QTimer::singleShot(0, [&] {
        auto* const dialog = qobject_cast<QInputDialog*>(QApplication::activeModalWidget());
        if (dialog == nullptr) {
            return;
        }
        label = dialog->labelText();
        value = dialog->intValue();
        maximum = dialog->intMaximum();
        dialog->reject();  // and no file dialog
    });
    child<QAction>(window, u"print_report").trigger();
    CHECK(label == u"30 Pairs found.  Print best:"_s);
    CHECK(value == 20);
    CHECK(maximum == 9999);
}
