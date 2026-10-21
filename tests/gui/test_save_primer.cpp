#include "core/io/pd2.hpp"
#include "gui/main_window.hpp"
#include "gui_primer_helpers.hpp"
#include "gui_test_helpers.hpp"

#include <QAction>
#include <QByteArray>
#include <QCheckBox>
#include <QFile>
#include <QIODevice>
#include <QLabel>
#include <QLineEdit>
#include <QString>
#include <QStringView>
#include <QTemporaryDir>
#include <QTest>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <span>

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

}  // namespace

// Reference: primers/M13EXT.PD2, written by PRIMER.exe's Save Active
// Primer (tests/legacy/README.md, "Primer files .PD2", recorded
// 2026-10-18): Direct Entry GTAAAACGACGGCCAGT, 5' extension GAATTC,
// built-in criteria
TEST_CASE("Primer/Create -> Save Active Primer writes PRIMER.exe's .PD2 file", "[gui][pd2][legacy]") {
    auto window = MainWindow{};
    direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    child<QAction>(window, u"modify_extension").trigger();
    auto& extension = child<QLineEdit>(window, u"extension_field");
    QTest::keyClicks(&extension, u"GAATTC"_s);
    QTest::keyClick(&extension, Qt::Key_Return);
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"M13EXT.PD2"_s);
    REQUIRE_FALSE(window.save_primer_to(path).has_value());
    CHECK(read_all(path) == read_all(legacy_file(u"primers/M13EXT.PD2")));
}

// Same reference: primers/COMP.PD2 (TEST600.GB, 400, complement strand);
// the molecule's path is this computer's, not PRIMER.exe's C:\TEST600.GB
TEST_CASE("Save Active Primer: a primer made from a molecule", "[gui][pd2][legacy]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"pairs/TEST600.GB")).has_value());
    child<QAction>(window, u"from_molecule").trigger();
    child<QCheckBox>(window, u"complement_box").setChecked(true);
    auto& position = child<QLineEdit>(window, u"position_field");
    QTest::keyClicks(&position, u"400"_s);
    QTest::keyClick(&position, Qt::Key_Return);
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"COMP.PD2"_s);
    REQUIRE_FALSE(window.save_primer_to(path).has_value());
    auto const saved = read_all(path);
    auto const parsed = primerpicker::io::parse_pd2(std::as_bytes(std::span{saved.constData(), static_cast<std::size_t>(saved.size())}));
    auto const primer = parsed.primer.value_or(primerpicker::io::Pd2Primer{});
    CHECK(primer.sequence == "TGTATCCAGGTAATGCAAGG");
    CHECK(primer.molecule_name == "TEST600");
    CHECK(primer.position == 400);
    CHECK(primer.complement);
    // and it can be retrieved
    REQUIRE_FALSE(window.load_primer(path).has_value());
    CHECK(child<QLabel>(window, u"sequence_label").text() == u"5'-TGTATCCAGGTAATGCAAGG-3'"_s);
}

TEST_CASE("Save Active Primer without a primer writes nothing", "[gui][pd2]") {
    auto window = MainWindow{};
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"NONE.PD2"_s);
    CHECK(window.save_primer_to(path).has_value());
    CHECK_FALSE(QFile::exists(path));
}
