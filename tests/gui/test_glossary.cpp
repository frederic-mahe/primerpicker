#include "core/glossary/glossary.hpp"
#include "core/glossary/glossary_table.hpp"
#include "gui/main_window.hpp"
#include "gui_primer_helpers.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QDialog>
#include <QKeySequence>
#include <QListView>
#include <QString>
#include <QTableView>
#include <QTextBrowser>
#include <QVariant>

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstddef>
#include <string_view>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;
using primerpicker::gui::test::direct_entry;

auto from_view(std::string_view const text) -> QString {
    return QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size()));
}

}  // namespace

TEST_CASE("Help -> Glossary lists the built-in definitions (F1)", "[gui][glossary]") {
    auto window = MainWindow{};
    auto& glossary = child<QAction>(window, u"glossary");
    CHECK(glossary.shortcut() == QKeySequence{Qt::Key_F1});
    glossary.trigger();
    CHECK(child<QDialog>(window, u"glossary_dialog").isVisible());
    auto const& list = *child<QListView>(window, u"glossary_list").model();
    REQUIRE(list.rowCount() == static_cast<int>(primerpicker::glossary::glossary_table.size()));
    CHECK(list.data(list.index(0, 0)).toString() == u"Primer"_s);
}

// Reference: data/glossary/glossary.md (reviewed by the lab, 2026-10-21)
TEST_CASE("the glossary shows an entry's text", "[gui][glossary]") {
    auto window = MainWindow{};
    window.show_glossary(u"Tm"_s);
    auto const& text = child<QTextBrowser>(window, u"glossary_text");
    CHECK(text.toPlainText().contains(u"The temperature at which half of the primer molecules"_s));
    CHECK(text.toPlainText().startsWith(u"Melting temperature (Tm) of a primer"_s));
}

TEST_CASE("the criteria table explains each criterion in a tooltip", "[gui][glossary]") {
    auto window = MainWindow{};
    direct_entry(window, u"GTAAAACGACGGCCAGT"_s);
    auto const& model = *child<QTableView>(window, u"criteria_table").model();
    constexpr auto criterion_column = 0;
    auto const tooltip = [&model](int const row) {
        return model.data(model.index(row, criterion_column), Qt::ToolTipRole).toString();
    };
    auto const definition = [](std::string_view const name) {
        auto const entry = primerpicker::glossary::find_entry(name);
        return entry.has_value() ? from_view(primerpicker::glossary::definition(*entry)) : QString{};
    };
    constexpr auto entries =
        std::array<std::string_view, 7>{"gc", "tm", "hairpin", "dimer", "dimer", "runs", "gc-clamp"};
    for (auto row = std::size_t{0}; row < entries.size(); ++row) {
        CHECK(tooltip(static_cast<int>(row)) == definition(entries.at(row)));
    }
}
