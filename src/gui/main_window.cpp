#include "gui/main_window.hpp"

#include "core/io/molecule_file.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "core/sequence/topology.hpp"
#include "gui/colour_scheme.hpp"

#include <QAction>
#include <QAnyStringView>
#include <QApplication>
#include <QDialog>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QIODevice>
#include <QMessageBox>
#include <QMenu>
#include <QMenuBar>
#include <QObject>
#include <QString>
#include <QStatusBar>
#include <QStringView>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

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

constexpr auto colour_scheme_key = QAnyStringView{u"appearance/colour_scheme"};

auto topology_name(sequence::Topology const topology) -> QString {
    switch (topology) {
        case sequence::Topology::linear:
            return u"linear"_s;
        case sequence::Topology::circular:
            return u"circular"_s;
        default:
            return {};  // unreachable: all enumerators are handled
    }
}

}  // namespace

MainWindow::MainWindow() : colour_schemes_{this}, direct_entry_{*this}, from_molecule_{*this} {
    setWindowTitle(u"PrimerPicker"_s);
    add_legacy_menus();
    no_primer_.setObjectName(u"no_primer_label"_s);
    no_primer_.setText(QString::fromUtf16(no_primer_message));
    no_primer_.setAlignment(Qt::AlignCenter);
    pages_.addWidget(&no_primer_);
    pages_.addWidget(&analysis_);
    setCentralWidget(&pages_);
    molecule_label_.setObjectName(u"molecule_label"_s);
    statusBar()->addWidget(&molecule_label_);
    settings_label_.setObjectName(u"settings_label"_s);
    statusBar()->addPermanentWidget(&settings_label_);
    QObject::connect(&direct_entry_, &QDialog::accepted, this, &MainWindow::use_direct_entry);
    QObject::connect(&from_molecule_, &FromMoleculeDialog::position_entered, this, &MainWindow::use_from_molecule);
}

MainWindow::MainWindow(QSettings& settings) : MainWindow() {
    settings_ = settings;
    auto const scheme = colour_scheme_from_name(settings.value(colour_scheme_key).toString());
    for (auto* const action : colour_schemes_.actions()) {
        action->setChecked(action->data().toString() == colour_scheme_name(scheme));
    }
    apply_colour_scheme(scheme);
}

// The menus of PRIMER.exe (proposal §1.1); & marks the Alt+letter
// mnemonics, "..." the items that open a dialog
void MainWindow::add_legacy_menus() {
    auto& menu_bar = *menuBar();
    auto& primer_create = *menu_bar.addMenu(u"&Primer/Create"_s);
    auto const& from_molecule = add_item(primer_create, {.text = u"From &Molecule...", .name = u"from_molecule"});
    QObject::connect(&from_molecule, &QAction::triggered, this, &MainWindow::create_from_molecule);
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
    // not in PRIMER.exe: the window's colours (proposal §3.2)
    modify.addSeparator();
    add_colour_scheme_menu(modify);

    auto& analyze = *menu_bar.addMenu(u"&Analyze"_s);
    auto const& active_primer = add_item(analyze, {.text = u"&Active Primer", .name = u"analyze_active_primer"});
    QObject::connect(&active_primer, &QAction::triggered, this, &MainWindow::analyze_active_primer);
    add_disabled_items(analyze, {u"Primer &Pair"_s, u"&Multi Primer"_s});

    auto& molecule_ops = *menu_bar.addMenu(u"M&olecule-Ops"_s);
    add_disabled_items(molecule_ops, {u"Find &Enzymes"_s, u"&View/Edit"_s, u"&Save/Remove"_s});

    auto const& quit = *menu_bar.addAction(u"&Quit"_s);
    QObject::connect(&quit, &QAction::triggered, &QApplication::closeAllWindows);
}

void MainWindow::add_colour_scheme_menu(QMenu& parent) {
    auto& menu = *parent.addMenu(u"Colour Sc&heme"_s);
    struct SchemeItem {
        QStringView text;
        ColourScheme scheme = ColourScheme::system;
    };
    constexpr auto items = std::array{SchemeItem{.text = u"&System", .scheme = ColourScheme::system},
                                      SchemeItem{.text = u"&Light", .scheme = ColourScheme::light},
                                      SchemeItem{.text = u"&Dark", .scheme = ColourScheme::dark}};
    for (auto const& item : items) {
        auto& action = *menu.addAction(item.text.toString());
        auto const name = colour_scheme_name(item.scheme);
        action.setObjectName(u"colour_scheme_"_s + name);
        action.setData(name);
        action.setCheckable(true);
        action.setChecked(item.scheme == ColourScheme::system);
        colour_schemes_.addAction(&action);
        QObject::connect(&action, &QAction::triggered, this, [this, scheme = item.scheme] {
            choose_colour_scheme(scheme);
        });
    }
}

void MainWindow::choose_colour_scheme(ColourScheme const scheme) {
    apply_colour_scheme(scheme);
    if (settings_.has_value()) {
        settings_->get().setValue(colour_scheme_key, colour_scheme_name(scheme));
    }
}

auto MainWindow::load_settings(QString const& path) -> std::optional<QString> {
    auto file = QFile{path};
    if (not file.open(QIODevice::ReadOnly)) {
        return u"cannot read %1"_s.arg(path);
    }
    auto const contents = file.readAll();
    auto bytes = std::vector<std::byte>{};
    bytes.reserve(static_cast<std::size_t>(contents.size()));
    std::ranges::transform(contents, std::back_inserter(bytes),
                           [](char const character) { return std::bit_cast<std::byte>(character); });
    auto const parsed = io::parse_pd2_cfg(bytes);
    if (not parsed.defaults.has_value()) {
        return u"%1: not a Primer Designer 2.0 settings file (pd2.cfg)"_s.arg(path);
    }
    defaults_ = *parsed.defaults;
    settings_label_.setText(u"Settings: %1"_s.arg(QFileInfo{path}.fileName()));
    return std::nullopt;
}

auto MainWindow::load_molecule(QString const& path) -> std::optional<QString> {
    auto file = QFile{path};
    if (not file.open(QIODevice::ReadOnly)) {
        return u"cannot read %1"_s.arg(path);
    }
    auto const bytes = file.readAll();
    auto const name = QFileInfo{path}.completeBaseName().toStdString();
    auto const parsed = io::parse_molecule_file(
        {.text = std::string_view{bytes.constData(), static_cast<std::size_t>(bytes.size())}, .name = name});
    auto sequence = parsed.molecule.has_value() ? sequence::to_nucleotides(parsed.molecule->sequence) : std::nullopt;
    if (not sequence.has_value()) {
        return u"%1: not a DNA molecule (GenBank, FASTA or plain sequence of A, C, G and T)"_s.arg(path);
    }
    molecule_ = LoadedMolecule{.name = QString::fromStdString(parsed.molecule->name),
                               .sequence = std::move(*sequence),
                               .topology = parsed.molecule->topology};
    molecule_label_.setText(u"Molecule: %1, %2 bp, %3"_s.arg(molecule_->name)
                                .arg(molecule_->sequence.size())
                                .arg(topology_name(molecule_->topology)));
    return std::nullopt;
}

// Load Molecule from Disk: a file chosen in a file dialog; whether a
// molecule was loaded
auto MainWindow::ask_for_molecule() -> bool {
    auto const path = QFileDialog::getOpenFileName(this, u"Load Molecule from Disk"_s, {},
                                                   u"Molecule files (*.gb *.gbk *.fa *.fasta *.seq *.txt);;"
                                                   u"All files (*)"_s);
    if (path.isEmpty()) {
        return false;
    }
    if (auto const error = load_molecule(path); error.has_value()) {
        QMessageBox::warning(this, u"Load Molecule"_s, *error);
        return false;
    }
    return true;
}

// As PRIMER.exe, a loaded molecule is used again without asking
void MainWindow::create_from_molecule() {
    if (not molecule_.has_value() and not ask_for_molecule()) {
        return;
    }
    if (molecule_.has_value()) {  // always, here
        from_molecule_.start(molecule_->name);
    }
}

void MainWindow::use_from_molecule() {
    auto const position = from_molecule_.position();
    auto primer = position.has_value() and molecule_.has_value()
                      ? sequence::primer_at(molecule_->sequence, {.position = *position,
                                                                  .length = defaults_.primer_length,
                                                                  .strand = from_molecule_.strand()})
                      : std::nullopt;
    if (not primer.has_value()) {
        from_molecule_.refuse_position();
        return;
    }
    from_molecule_.accept();
    active_primer_ = ActivePrimer{.sequence = std::move(*primer), .extension_length = 0};
    analyze_active_primer();
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
