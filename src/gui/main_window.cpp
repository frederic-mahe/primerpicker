#include "gui/main_window.hpp"

#include <QAction>
#include <QApplication>
#include <QMenu>
#include <QMenuBar>
#include <QObject>
#include <QString>

#include <initializer_list>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// Disabled items: features of Primer Designer not implemented in the
// GUI yet
void add_disabled_items(QMenu& menu, std::initializer_list<QString> const texts) {
    for (auto const& text : texts) {
        menu.addAction(text)->setEnabled(false);
    }
}

// The menus of PRIMER.exe (proposal §1.1); & marks the Alt+letter
// mnemonics, "..." the items that open a dialog
void add_legacy_menus(QMenuBar& menu_bar) {
    auto& primer_create = *menu_bar.addMenu(u"&Primer/Create"_s);
    add_disabled_items(primer_create, {u"From &Molecule..."_s, u"&Direct Entry..."_s, u"Match &Active Primer"_s,
                                       u"Find &Pair"_s, u"&Save Active Primer"_s, u"&Retrieve from disk"_s});

    auto& modify = *menu_bar.addMenu(u"&Modify"_s);
    modify.addSection(u"Active Primer:"_s);
    add_disabled_items(modify,
                       {u"&Position"_s, u"&Length"_s, u"5'-&extension"_s, u"&Criteria"_s, u"M&utagenize"_s});
    modify.addSeparator();
    add_disabled_items(modify, {u"&System Defaults"_s});

    auto& analyze = *menu_bar.addMenu(u"&Analyze"_s);
    add_disabled_items(analyze, {u"&Active Primer"_s, u"Primer &Pair"_s, u"&Multi Primer"_s});

    auto& molecule_ops = *menu_bar.addMenu(u"M&olecule-Ops"_s);
    add_disabled_items(molecule_ops, {u"Find &Enzymes"_s, u"&View/Edit"_s, u"&Save/Remove"_s});

    auto const& quit = *menu_bar.addAction(u"&Quit"_s);
    QObject::connect(&quit, &QAction::triggered, &QApplication::closeAllWindows);
}

}  // namespace

MainWindow::MainWindow() {
    setWindowTitle(u"PrimerPicker"_s);
    add_legacy_menus(*menuBar());
}

}  // namespace primerpicker::gui
