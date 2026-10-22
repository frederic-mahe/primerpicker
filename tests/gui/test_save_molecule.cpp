#include "gui/main_window.hpp"

#include "gui_test_helpers.hpp"

#include <QByteArray>
#include <QFile>
#include <QIODevice>
#include <QLabel>
#include <QString>
#include <QStringView>
#include <QTemporaryDir>

#include <catch2/catch_test_macros.hpp>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

auto bytes_of(QString const& path) -> QByteArray {
    auto file = QFile{path};
    REQUIRE(file.open(QIODevice::ReadOnly));
    return file.readAll();
}

}  // namespace

// Reference: PRIMER.exe, Molecule-Ops -> Save/Remove -> S
// (tests/legacy/molecule_save/, recorded 2026-10-22): TEST600.GB saved
// as ZZ.GB, TCIRC.GB (circular) as ZC.GB, TFASTA.FA over itself, in
// PRIMER.exe's own format
TEST_CASE("Molecule-Ops -> Save/Remove: Save writes PRIMER.exe's format", "[gui][molecule]") {
    struct Case {
        QStringView loaded;
        QStringView written;
    };
    for (auto const& [loaded, written] :
         {Case{u"pairs/TEST600.GB", u"molecule_save/written/save_new.ZZ.GB"},
          Case{u"circular/TCIRC.GB", u"molecule_save/written/save_circ.ZC.GB"},
          Case{u"molecules/TFASTA.FA", u"molecule_save/written/save_fa.TFASTA.FA"}}) {
        CAPTURE(loaded.toString());
        auto window = MainWindow{};
        REQUIRE_FALSE(window.load_molecule(legacy_file(loaded)).has_value());
        auto const directory = QTemporaryDir{};
        REQUIRE(directory.isValid());
        auto const path = directory.filePath(u"SAVED"_s);
        REQUIRE_FALSE(window.save_molecule_to(path).has_value());
        CHECK(bytes_of(path) == bytes_of(legacy_file(written)));
    }
}

// Reference: PRIMER.exe, "Molecule List is Empty" (view_edit/
// saveremove_open, recorded 2026-10-21)
TEST_CASE("Molecule-Ops -> Save/Remove: Save needs a molecule", "[gui][molecule]") {
    auto window = MainWindow{};
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    CHECK(window.save_molecule_to(directory.filePath(u"SAVED"_s)) == u"Molecule List is Empty"_s);
}

// Reference: PRIMER.exe loads the files it writes (molecule_save/
// reload_header: ZZ.GB is TEST600, 600 bases) and refuses malformed ones
// with "Illegal File Format" (molecule_save/hand9)
TEST_CASE("PRIMER.exe's own molecule files load", "[gui][molecule]") {
    auto window = MainWindow{};
    REQUIRE_FALSE(window.load_molecule(legacy_file(u"molecule_save/written/save_new.ZZ.GB")).has_value());
    REQUIRE_FALSE(window.view_molecule().has_value());
    CHECK(child<QLabel>(window, u"molecule_name_label").text() == u"TEST600"_s);
    auto const error = window.load_molecule(legacy_file(u"molecule_save/inputs/ZH9"));
    CHECK(error.value_or(QString{}).contains(u"Illegal File Format"_s));
}
