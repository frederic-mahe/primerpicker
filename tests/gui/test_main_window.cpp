#include "gui/main_window.hpp"

#include <QAction>
#include <QList>
#include <QMenu>
#include <QMenuBar>
#include <QString>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;

auto texts(QList<QAction*> const& actions) -> QList<QString> {
    auto result = QList<QString>{};
    for (auto const* const action : actions) {
        result.append(action->isSeparator() ? u"---"_s : action->text());
    }
    return result;
}

auto menu_titled(MainWindow const& window, QString const& title) -> QMenu const& {
    auto const actions = window.menuBar()->actions();
    auto const found = std::ranges::find_if(actions, [&title](QAction const* const action) {
        return action->text() == title and action->menu() != nullptr;
    });
    REQUIRE(found != actions.end());
    return *(*found)->menu();
}

}  // namespace

// Reference: the main menu of PRIMER.exe (proposal §1.1, from the
// strings of examples/PRIMER.exe); & marks the Alt+letter mnemonics
TEST_CASE("the menu bar mirrors Primer Designer's main menu, with a View menu", "[gui][menu]") {
    auto const window = MainWindow{};
    CHECK(texts(window.menuBar()->actions()) ==
          QList<QString>{u"&Primer/Create"_s, u"&Modify"_s, u"&Analyze"_s, u"M&olecule-Ops"_s, u"&View"_s, u"&Help"_s,
                         u"&Quit"_s});
}

TEST_CASE("each menu lists Primer Designer's items in order", "[gui][menu]") {
    auto const window = MainWindow{};
    CHECK(texts(menu_titled(window, u"&Primer/Create"_s).actions()) ==
          QList<QString>{u"From &Molecule..."_s, u"&Direct Entry..."_s, u"Match &Active Primer..."_s, u"Find &Pair..."_s,
                      u"&Save Active Primer..."_s, u"Save to &List"_s, u"&Retrieve from disk..."_s});
    // the section title "Active Primer:" is a separator with a text
    CHECK(texts(menu_titled(window, u"&Modify"_s).actions()) ==
          QList<QString>{u"---"_s, u"&Position"_s, u"&Length"_s, u"5'-&extension..."_s, u"&Criteria..."_s,
                      u"M&utagenize"_s, u"---"_s, u"&System Defaults..."_s});
    CHECK(texts(menu_titled(window, u"&Analyze"_s).actions()) ==
          QList<QString>{u"&Active Primer"_s, u"Primer &Pair..."_s, u"&Multi Primer..."_s, u"---"_s,
                         u"Print to &Disk..."_s});
    CHECK(texts(menu_titled(window, u"M&olecule-Ops"_s).actions()) ==
          QList<QString>{u"Find &Enzymes..."_s, u"&View/Edit"_s, u"&Save/Remove"_s});
    // not in PRIMER.exe: the window's own settings
    CHECK(texts(menu_titled(window, u"&View"_s).actions()) ==
          QList<QString>{u"Primer &List"_s, u"&Colour Scheme"_s});
}

TEST_CASE("menu items that are not implemented yet are disabled", "[gui][menu]") {
    auto const window = MainWindow{};
    auto enabled = QList<QString>{};
    for (auto const* const menu_action : window.menuBar()->actions()) {
        auto const* const menu = menu_action->menu();
        if (menu == nullptr) {
            continue;
        }
        for (auto const* const action : menu->actions()) {
            if (not action->isSeparator() and action->isEnabled()) {
                enabled.append(action->text());
            }
        }
    }
    CHECK(enabled == QList<QString>{u"From &Molecule..."_s, u"&Direct Entry..."_s, u"Match &Active Primer..."_s,
                                    u"Find &Pair..."_s, u"&Save Active Primer..."_s,
                                    u"Save to &List"_s, u"&Retrieve from disk..."_s,
                                    u"5'-&extension..."_s, u"&Criteria..."_s, u"&System Defaults..."_s,
                                    u"&Active Primer"_s, u"Primer &Pair..."_s, u"&Multi Primer..."_s,
                                    u"Print to &Disk..."_s,
                                    u"Find &Enzymes..."_s, u"Primer &List"_s, u"&Colour Scheme"_s,
                                    u"&Glossary"_s});
}
