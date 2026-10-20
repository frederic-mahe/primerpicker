#include "gui/main_window.hpp"

#include <QAction>
#include <QApplication>
#include <QDialog>
#include <QMenu>
#include <QMenuBar>
#include <QObject>
#include <QString>
#include <QStringView>

#include <initializer_list>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// PRIMER.exe's message when there is no active primer (from the strings
// of examples/PRIMER.exe)
constexpr auto no_primer_message = u"Please create or retrieve a primer:\n\n"
                                   u"Select From Molecule to create a primer using a sequence file.\n"
                                   u"Select Direct Entry to type in the sequence of a new primer.";

// An item of a menu: its text and its objectName. Views passed by
// value, not QString temporaries: GCC warns that a reference returned
// by add_item might dangle.
struct MenuItem {
    QStringView text;
    QStringView name;
};

auto add_item(QMenu& menu, MenuItem const item) -> QAction& {
    auto& action = *menu.addAction(item.text.toString());
    action.setObjectName(item.name);
    return action;
}

// Disabled items: features of Primer Designer not implemented in the
// GUI yet
void add_disabled_items(QMenu& menu, std::initializer_list<QString> const texts) {
    for (auto const& text : texts) {
        menu.addAction(text)->setEnabled(false);
    }
}

}  // namespace

MainWindow::MainWindow() : direct_entry_{*this} {
    setWindowTitle(u"PrimerPicker"_s);
    add_legacy_menus();
    no_primer_.setObjectName(u"no_primer_label"_s);
    no_primer_.setText(QString::fromUtf16(no_primer_message));
    no_primer_.setAlignment(Qt::AlignCenter);
    pages_.addWidget(&no_primer_);
    pages_.addWidget(&analysis_);
    setCentralWidget(&pages_);
    QObject::connect(&direct_entry_, &QDialog::accepted, this, &MainWindow::use_direct_entry);
}

// The menus of PRIMER.exe (proposal §1.1); & marks the Alt+letter
// mnemonics, "..." the items that open a dialog
void MainWindow::add_legacy_menus() {
    auto& menu_bar = *menuBar();
    auto& primer_create = *menu_bar.addMenu(u"&Primer/Create"_s);
    add_disabled_items(primer_create, {u"From &Molecule..."_s});
    auto const& direct_entry = add_item(primer_create, {.text = u"&Direct Entry...", .name = u"direct_entry"});
    QObject::connect(&direct_entry, &QAction::triggered, &direct_entry_, &DirectEntryDialog::start);
    add_disabled_items(primer_create, {u"Match &Active Primer"_s, u"Find &Pair"_s, u"&Save Active Primer"_s,
                                       u"&Retrieve from disk"_s});

    auto& modify = *menu_bar.addMenu(u"&Modify"_s);
    modify.addSection(u"Active Primer:"_s);
    add_disabled_items(modify,
                       {u"&Position"_s, u"&Length"_s, u"5'-&extension"_s, u"&Criteria"_s, u"M&utagenize"_s});
    modify.addSeparator();
    add_disabled_items(modify, {u"&System Defaults"_s});

    auto& analyze = *menu_bar.addMenu(u"&Analyze"_s);
    auto const& active_primer = add_item(analyze, {.text = u"&Active Primer", .name = u"analyze_active_primer"});
    QObject::connect(&active_primer, &QAction::triggered, this, &MainWindow::analyze_active_primer);
    add_disabled_items(analyze, {u"Primer &Pair"_s, u"&Multi Primer"_s});

    auto& molecule_ops = *menu_bar.addMenu(u"M&olecule-Ops"_s);
    add_disabled_items(molecule_ops, {u"Find &Enzymes"_s, u"&View/Edit"_s, u"&Save/Remove"_s});

    auto const& quit = *menu_bar.addAction(u"&Quit"_s);
    QObject::connect(&quit, &QAction::triggered, &QApplication::closeAllWindows);
}

void MainWindow::use_direct_entry() {
    active_primer_ = ActivePrimer{.sequence = direct_entry_.primer(), .extension_length = 0};
    analyze_active_primer();
}

void MainWindow::analyze_active_primer() {
    if (not active_primer_.has_value()) {
        pages_.setCurrentWidget(&no_primer_);
        return;
    }
    analysis_.show_primer({.sequence = active_primer_->sequence, .extension_length = active_primer_->extension_length},
                          defaults_.criteria);
    pages_.setCurrentWidget(&analysis_);
}

}  // namespace primerpicker::gui
