#include "gui/main_window.hpp"

#include "core/analysis/legacy_homology.hpp"
#include "core/analysis/legacy_multi_primer.hpp"
#include "core/enzymes/rebase_enzymes.hpp"
#include "core/enzymes/site_search.hpp"
#include "core/io/legacy_native_molecule.hpp"
#include "core/io/molecule_file.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/search/legacy_match.hpp"
#include "core/search/legacy_picking.hpp"
#include "core/io/pd2.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/report/legacy_report.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "core/sequence/topology.hpp"
#include "gui/analysis_view.hpp"
#include "gui/colour_scheme.hpp"
#include "gui/find_enzymes_dialog.hpp"
#include "gui/find_pair_dialog.hpp"
#include "gui/homology_view.hpp"
#include "gui/match_dialog.hpp"
#include "gui/molecule_view.hpp"
#include "gui/multi_primer_dialog.hpp"
#include "gui/mutagenize_view.hpp"
#include "gui/pair_dialog.hpp"
#include "gui/picking_view.hpp"
#include "gui/sequence_dialog.hpp"

#include <QAbstractButton>
#include <QAbstractItemView>
#include <QAction>
#include <QAnyStringView>
#include <QChar>
#include <QDate>
#include <QApplication>
#include <QByteArray>
#include <QDialog>
#include <QFutureWatcherBase>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QInputDialog>
#include <QFileInfo>
#include <QHeaderView>
#include <QIODevice>
#include <QKeySequence>
#include <QList>
#include <QMessageBox>
#include <QModelIndex>
#include <QProgressDialog>
#include <QPushButton>
#include <QMenu>
#include <QMenuBar>
#include <QObject>
#include <QString>
#include <QStatusBar>
#include <QStringView>
#include <QtConcurrentRun>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstddef>
#include <iterator>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
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

constexpr auto colour_scheme_key = QAnyStringView{u"appearance/colour_scheme"};

auto letters(std::span<sequence::Nucleotide const> const nucleotides) -> QString {
    auto text = QString{};
    for (auto const nucleotide : nucleotides) {
        text.append(QChar::fromLatin1(sequence::to_letter(nucleotide)));
    }
    return text;
}

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

// the date of the reports
auto today() -> std::chrono::year_month_day {
    auto const date = QDate::currentDate();
    return std::chrono::year_month_day{std::chrono::year{date.year()},
                                       std::chrono::month{static_cast<unsigned>(date.month())},
                                       std::chrono::day{static_cast<unsigned>(date.day())}};
}

// the error, if REPORT could not be written to PATH
auto write_report(QString const& path, std::string const& report) -> std::optional<QString> {
    auto file = QFile{path};
    if (not file.open(QIODevice::WriteOnly) or
        file.write(report.data(), static_cast<qint64>(report.size())) != static_cast<qint64>(report.size())) {
        return u"cannot write the report to %1"_s.arg(path);
    }
    return std::nullopt;
}

// a mutated primer: PRIMER, edited on the Mutagenize screen
auto mutated_copy(ActivePrimer const& primer, std::vector<sequence::Nucleotide> mutated) -> ActivePrimer {
    auto copy = primer;
    copy.original = primer.sequence;
    copy.sequence = std::move(mutated);
    copy.origin = primer.origin + u" (mutated)"_s;
    return copy;
}

}  // namespace

MainWindow::MainWindow()
    : list_model_{{u"#"_s, u"Sequence"_s, u"Length"_s, u"Origin"_s}},
      colour_schemes_{this},
      direct_entry_{*this, direct_entry_kind},
      extension_{*this, extension_kind},
      from_molecule_{*this},
      find_enzymes_{*this},
      find_pair_{*this},
      pair_dialog_{*this},
      multi_dialog_{*this},
      match_dialog_{*this},
      find_pair_progress_{this},
      glossary_{*this},
      criteria_dialog_{*this},
      defaults_dialog_{*this} {
    setWindowTitle(u"PrimerPicker"_s);
    // the dock's title first: it names the dock's View menu item
    list_dock_.setWindowTitle(u"Primer List"_s);
    add_legacy_menus();
    no_primer_.setObjectName(u"no_primer_label"_s);
    no_primer_.setText(QString::fromUtf16(no_primer_message));
    no_primer_.setAlignment(Qt::AlignCenter);
    pages_.addWidget(&no_primer_);
    pages_.addWidget(&analysis_);
    pages_.addWidget(&sites_);
    pages_.addWidget(&pairs_);
    pages_.addWidget(&homology_);
    pages_.addWidget(&pair_);
    pages_.addWidget(&multi_);
    pages_.addWidget(&match_);
    pages_.addWidget(&picking_);
    pages_.addWidget(&mutagenize_);
    pages_.addWidget(&molecule_view_);
    list_view_.setObjectName(u"primer_list"_s);
    list_view_.setModel(&list_model_);
    list_view_.setEditTriggers(QAbstractItemView::NoEditTriggers);
    list_view_.setSelectionBehavior(QAbstractItemView::SelectRows);
    list_view_.setSelectionMode(QAbstractItemView::SingleSelection);
    list_view_.verticalHeader()->hide();
    list_view_.horizontalHeader()->setStretchLastSection(true);
    list_dock_.setObjectName(u"primer_list_dock"_s);
    list_dock_.setWidget(&list_view_);
    addDockWidget(Qt::LeftDockWidgetArea, &list_dock_);
    QObject::connect(&list_view_, &QAbstractItemView::activated, this,
                     [this](QModelIndex const& index) { show_listed_primer(index.row()); });
    QObject::connect(&pair_dialog_, &QDialog::accepted, this, &MainWindow::use_primer_pair);
    QObject::connect(&multi_dialog_, &QDialog::accepted, this, &MainWindow::use_multi_primer);
    QObject::connect(&match_dialog_, &MatchDialog::search_entered, this, &MainWindow::use_match);
    QObject::connect(&criteria_dialog_, &CriteriaDialog::criteria_entered, this, &MainWindow::use_criteria);
    QObject::connect(&defaults_dialog_, &SystemDefaultsDialog::defaults_entered, this,
                     &MainWindow::use_system_defaults);
    QObject::connect(&defaults_dialog_, &SystemDefaultsDialog::load_requested, this, &MainWindow::load_settings_file);
    QObject::connect(&defaults_dialog_, &SystemDefaultsDialog::save_requested, this, &MainWindow::save_settings_file);
    setCentralWidget(&pages_);
    molecule_label_.setObjectName(u"molecule_label"_s);
    statusBar()->addWidget(&molecule_label_);
    settings_label_.setObjectName(u"settings_label"_s);
    statusBar()->addPermanentWidget(&settings_label_);
    QObject::connect(&direct_entry_, &QDialog::accepted, this, &MainWindow::use_direct_entry);
    QObject::connect(&extension_, &QDialog::accepted, this, &MainWindow::use_extension);
    QObject::connect(&find_enzymes_, &FindEnzymesDialog::enzyme_entered, this, &MainWindow::use_find_enzymes);
    QObject::connect(&find_pair_, &FindPairDialog::search_entered, this, &MainWindow::use_find_pair);
    find_pair_progress_.setObjectName(u"find_pair_progress"_s);
    find_pair_progress_.setWindowTitle(u"Find Pair"_s);
    find_pair_progress_.setLabelText(u"Searching for primer pairs..."_s);
    find_pair_progress_.setRange(0, 0);  // busy: the search gives no progress
    find_pair_progress_.setWindowModality(Qt::WindowModal);
    find_pair_progress_.setAutoClose(false);
    find_pair_progress_.setAutoReset(false);
    find_pair_progress_.hide();
    QObject::connect(&find_pair_progress_, &QProgressDialog::canceled, this, &MainWindow::cancel_find_pair);
    QObject::connect(&find_pair_watcher_, &QFutureWatcherBase::finished, this, &MainWindow::find_pair_finished);
    QObject::connect(&analysis_, &AnalysisView::print_requested, this, &MainWindow::print_report);
    QObject::connect(&analysis_, &AnalysisView::uniqueness_requested, this,
                     [this] { show_homology(analysis::HomologyKind::uniqueness); });
    QObject::connect(&analysis_, &AnalysisView::false_priming_requested, this,
                     [this] { show_homology(analysis::HomologyKind::false_priming); });
    QObject::connect(&from_molecule_, &FromMoleculeDialog::position_entered, this, &MainWindow::use_from_molecule);
    QObject::connect(&picking_, &PickingView::location_changed, this, &MainWindow::use_picked_location);
    QObject::connect(&mutagenize_, &MutagenizeView::save_requested, this, &MainWindow::save_mutated_primer);
    QObject::connect(&mutagenize_, &MutagenizeView::inspect_requested, this, &MainWindow::inspect_mutated_primer);
    QObject::connect(&mutagenize_, &MutagenizeView::exit_requested, this, &MainWindow::analyze_active_primer);
    // S, Enter or Esc: as PRIMER.exe, the primer joins the Primer List and
    // the picking screen closes
    QObject::connect(&picking_, &PickingView::save_requested, this, &MainWindow::save_to_list);
}

// a running search reads the molecule: it is stopped first
MainWindow::~MainWindow() {
    cancel_find_pair();
    find_pair_watcher_.waitForFinished();
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
    QObject::connect(&direct_entry, &QAction::triggered, &direct_entry_, [this] { direct_entry_.start(); });
    auto const& match = add_item(primer_create, {.text = u"Match &Active Primer...", .name = u"match_active_primer"});
    QObject::connect(&match, &QAction::triggered, this, &MainWindow::open_match);
    auto const& find_pair_item = add_item(primer_create, {.text = u"Find &Pair...", .name = u"find_pair"});
    QObject::connect(&find_pair_item, &QAction::triggered, this, &MainWindow::open_find_pair);
    auto const& save_item = add_item(primer_create, {.text = u"&Save Active Primer...", .name = u"save_active_primer"});
    QObject::connect(&save_item, &QAction::triggered, this, &MainWindow::save_primer);
    // not in PRIMER.exe as a menu item: its S (Save+Exit) key
    auto& save_to_list_item = add_item(primer_create, {.text = u"Save to &List", .name = u"save_to_list"});
    save_to_list_item.setShortcut(QKeySequence{Qt::CTRL | Qt::Key_L});
    QObject::connect(&save_to_list_item, &QAction::triggered, this, &MainWindow::save_to_list);
    auto const& retrieve = add_item(primer_create, {.text = u"&Retrieve from disk...", .name = u"retrieve_primer"});
    QObject::connect(&retrieve, &QAction::triggered, this, &MainWindow::retrieve_primer);

    auto& modify = *menu_bar.addMenu(u"&Modify"_s);
    modify.addSection(u"Active Primer:"_s);
    auto const& position = add_item(modify, {.text = u"&Position", .name = u"modify_position"});
    QObject::connect(&position, &QAction::triggered, this, &MainWindow::open_modify_position);
    auto const& length = add_item(modify, {.text = u"&Length", .name = u"modify_length"});
    QObject::connect(&length, &QAction::triggered, this, &MainWindow::open_change_length);
    auto const& extension = add_item(modify, {.text = u"5'-&extension...", .name = u"modify_extension"});
    QObject::connect(&extension, &QAction::triggered, this, &MainWindow::modify_extension);
    auto const& criteria = add_item(modify, {.text = u"&Criteria...", .name = u"modify_criteria"});
    QObject::connect(&criteria, &QAction::triggered, this, &MainWindow::modify_criteria);
    auto const& mutagenize_item = add_item(modify, {.text = u"M&utagenize", .name = u"modify_mutagenize"});
    QObject::connect(&mutagenize_item, &QAction::triggered, this, &MainWindow::open_mutagenize);
    modify.addSeparator();
    auto const& system_defaults = add_item(modify, {.text = u"&System Defaults...", .name = u"system_defaults"});
    QObject::connect(&system_defaults, &QAction::triggered, this, &MainWindow::open_system_defaults);

    auto& analyze = *menu_bar.addMenu(u"&Analyze"_s);
    auto const& active_primer = add_item(analyze, {.text = u"&Active Primer", .name = u"analyze_active_primer"});
    QObject::connect(&active_primer, &QAction::triggered, this, &MainWindow::analyze_active_primer);
    auto const& primer_pair = add_item(analyze, {.text = u"Primer &Pair...", .name = u"analyze_primer_pair"});
    QObject::connect(&primer_pair, &QAction::triggered, this, &MainWindow::open_primer_pair);
    auto const& multi_primer = add_item(analyze, {.text = u"&Multi Primer...", .name = u"analyze_multi_primer"});
    QObject::connect(&multi_primer, &QAction::triggered, this, &MainWindow::open_multi_primer);
    // PRIMER.exe's P on each screen
    analyze.addSeparator();
    auto& print = add_item(analyze, {.text = u"Print to &Disk...", .name = u"print_report"});
    print.setShortcut(QKeySequence{Qt::CTRL | Qt::Key_P});
    QObject::connect(&print, &QAction::triggered, this, &MainWindow::print_report);

    auto& molecule_ops = *menu_bar.addMenu(u"M&olecule-Ops"_s);
    auto const& find_enzymes = add_item(molecule_ops, {.text = u"Find &Enzymes...", .name = u"find_enzymes"});
    QObject::connect(&find_enzymes, &QAction::triggered, this, &MainWindow::open_find_enzymes);
    auto const& view_edit = add_item(molecule_ops, {.text = u"&View/Edit", .name = u"view_molecule"});
    QObject::connect(&view_edit, &QAction::triggered, this, &MainWindow::open_view_molecule);
    auto const& save_remove = add_item(molecule_ops, {.text = u"&Save/Remove", .name = u"save_remove_molecule"});
    QObject::connect(&save_remove, &QAction::triggered, this, &MainWindow::open_save_remove);

    // not in PRIMER.exe: the window's own settings, such as its colours
    // (proposal §3.2)
    auto& view = *menu_bar.addMenu(u"&View"_s);
    auto& list_toggle = *list_dock_.toggleViewAction();
    list_toggle.setText(u"Primer &List"_s);
    list_toggle.setObjectName(u"view_primer_list"_s);
    view.addAction(&list_toggle);
    add_colour_scheme_menu(view);

    // not in PRIMER.exe: the built-in definitions
    auto& help = *menu_bar.addMenu(u"&Help"_s);
    auto& glossary = add_item(help, {.text = u"&Glossary", .name = u"glossary"});
    glossary.setShortcut(QKeySequence{Qt::Key_F1});
    QObject::connect(&glossary, &QAction::triggered, this, [this] { show_glossary({}); });

    auto const& quit = *menu_bar.addAction(u"&Quit"_s);
    QObject::connect(&quit, &QAction::triggered, &QApplication::closeAllWindows);
}

void MainWindow::add_colour_scheme_menu(QMenu& parent) {
    auto& menu = *parent.addMenu(u"&Colour Scheme"_s);
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
    settings_name_ = QFileInfo{path}.fileName();
    settings_path_ = path;
    saved_defaults_ = defaults_;
    settings_label_.setText(u"Settings: %1"_s.arg(settings_name_));
    return std::nullopt;
}

auto MainWindow::save_settings_to(QString const& path) const -> std::optional<QString> {
    auto contents = QByteArray{};
    for (auto const byte : io::write_pd2_cfg(defaults_)) {
        contents.append(std::bit_cast<char>(byte));
    }
    auto file = QFile{path};
    if (not file.open(QIODevice::WriteOnly) or file.write(contents) != contents.size()) {
        return u"cannot write the settings to %1"_s.arg(path);
    }
    return std::nullopt;
}

void MainWindow::modify_criteria() {
    if (not active_primer_.has_value()) {
        analyze_active_primer();  // asks for a primer
        return;
    }
    criteria_dialog_.start(active_primer_->criteria);
}

// the active primer's own criteria (a .PD2 file keeps them)
void MainWindow::use_criteria() {
    if (active_primer_.has_value()) {
        active_primer_->criteria = criteria_dialog_.criteria();
        analyze_active_primer();
    }
}

void MainWindow::open_system_defaults() {
    defaults_dialog_.start(defaults_);
}

// new primers and searches use them; the active primer keeps its
// criteria (as a primer of PRIMER.exe)
void MainWindow::use_system_defaults() {
    defaults_ = defaults_dialog_.defaults();
    auto const name = settings_name_.isEmpty() ? u"built-in"_s : settings_name_;
    settings_label_.setText(u"Settings: %1 (changed)"_s.arg(name));
}

void MainWindow::load_settings_file() {
    auto const path = QFileDialog::getOpenFileName(this, u"Load Settings File"_s, {},
                                                   u"Settings files (*.cfg *.CFG);;All files (*)"_s);
    if (path.isEmpty()) {
        return;
    }
    if (auto const error = load_settings(path); error.has_value()) {
        QMessageBox::warning(this, u"Load Settings File"_s, *error);
        return;
    }
    defaults_dialog_.start(defaults_);
}

// the values of the dialog, as PRIMER.exe saves PD2.CFG
void MainWindow::save_settings_file() {
    auto const path = QFileDialog::getSaveFileName(this, u"Save Settings File"_s, u"PD2.CFG"_s,
                                                   u"Settings files (*.cfg *.CFG);;All files (*)"_s);
    if (path.isEmpty()) {
        return;
    }
    auto contents = QByteArray{};
    for (auto const byte : io::write_pd2_cfg(defaults_dialog_.defaults())) {
        contents.append(std::bit_cast<char>(byte));
    }
    auto file = QFile{path};
    if (not file.open(QIODevice::WriteOnly) or file.write(contents) != contents.size()) {
        QMessageBox::warning(this, u"Save Settings File"_s, u"cannot write the settings to %1"_s.arg(path));
        return;
    }
    saved_defaults_ = defaults_dialog_.defaults();
}

auto MainWindow::settings_changed() const -> bool {
    // the two settings as their files: what quitting would save
    return io::write_pd2_cfg(defaults_) != io::write_pd2_cfg(saved_defaults_);
}

auto MainWindow::settings_save_path() const -> QString {
    auto const file_name = u"PD2.CFG"_s;  // the only name PRIMER.exe loads
    if (settings_path_.isEmpty()) {
        return QDir::current().filePath(file_name);
    }
    return QFileInfo{settings_path_}.absoluteDir().filePath(file_name);
}

void MainWindow::closeEvent(QCloseEvent* const event) {
    if (not settings_changed()) {
        QMainWindow::closeEvent(event);
        return;
    }
    auto const answer = QMessageBox::question(
        this, u"Quit"_s, u"System Defaults Changed But Not Saved!\nSave current defaults?"_s,
        QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel, QMessageBox::Yes);
    if (answer == QMessageBox::No) {
        QMainWindow::closeEvent(event);
        return;
    }
    if (answer == QMessageBox::Yes) {
        auto const path = QFileDialog::getSaveFileName(this, u"Save Settings File"_s, settings_save_path(),
                                                       u"Settings files (*.cfg *.CFG);;All files (*)"_s);
        if (not path.isEmpty()) {
            if (auto const error = save_settings_to(path); error.has_value()) {
                QMessageBox::warning(this, u"Save Settings File"_s, *error);
            } else {
                saved_defaults_ = defaults_;
                QMainWindow::closeEvent(event);
                return;
            }
        }
    }
    event->ignore();  // Cancel, or no file saved: the window stays
}

auto MainWindow::load_molecule(QString const& path) -> std::optional<QString> {
    if (find_pair_running()) {
        return u"A Find Pair search is running"_s;  // it reads the molecule
    }
    auto file = QFile{path};
    if (not file.open(QIODevice::ReadOnly)) {
        return u"cannot read %1"_s.arg(path);
    }
    auto const bytes = file.readAll();
    auto const name = QFileInfo{path}.completeBaseName().toStdString();
    auto const parsed = io::parse_molecule_file(
        {.text = std::string_view{bytes.constData(), static_cast<std::size_t>(bytes.size())}, .name = name});
    if (parsed.error == io::MoleculeError::illegal_format) {
        return u"%1: Illegal File Format"_s.arg(path);  // PRIMER.exe's message
    }
    auto sequence = parsed.molecule.has_value() ? sequence::to_nucleotides(parsed.molecule->sequence) : std::nullopt;
    if (not sequence.has_value()) {
        return u"%1: not a DNA molecule (GenBank, FASTA or plain sequence of A, C, G and T)"_s.arg(path);
    }
    molecule_ = LoadedMolecule{.name = QString::fromStdString(parsed.molecule->name),
                               .path = QFileInfo{path}.absoluteFilePath(),
                               .text = bytes.toStdString(),
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

auto MainWindow::find_enzyme(QString const& name, enzymes::Mismatches const mismatches) -> std::optional<QString> {
    if (not molecule_.has_value()) {
        return u"Please load a molecule"_s;
    }
    auto const resolved = enzymes::resolve_enzyme(name.toStdString());
    if (not resolved.has_value()) {
        return u"Unknown enzyme: %1"_s.arg(name);
    }
    auto const pattern = enzymes::parse_site(resolved->enzyme.site);
    if (not pattern.has_value()) {
        return u"%1: unexpected recognition sequence"_s.arg(name);
    }
    auto const sites = enzymes::find_sites(molecule_->sequence, *pattern, mismatches, molecule_->topology);
    sites_.show_sites({.enzyme = resolved->enzyme, .mismatches = mismatches, .molecule = molecule_->sequence}, sites);
    report_.reset();  // not recorded
    pages_.setCurrentWidget(&sites_);
    return std::nullopt;
}

// The settings of a Find Pair search on the loaded molecule, or the
// message of PRIMER.exe's form for targets that do not fit
auto MainWindow::find_pair_settings(FindPairRequest const& request) const
    -> std::variant<search::FindPairSettings, QString> {
    if (not molecule_.has_value()) {
        return u"Please load a molecule"_s;
    }
    auto const settings = search::FindPairSettings{.target_a = request.target_a,
                                                   .start_offset_a = request.start_offset_a,
                                                   .target_b = request.target_b,
                                                   .end_offset_b = request.end_offset_b,
                                                   .length = defaults_.primer_length,
                                                   .gc_tolerance_percent = defaults_.gc_tolerance_percent,
                                                   .tm_tolerance_celsius = defaults_.tm_tolerance_celsius};
    // as PRIMER.exe's form, targets whose primers do not fit are refused
    if (auto const refused = search::legacy_check_targets(settings, molecule_->sequence.size());
        refused != search::FindPairTargetError::none) {
        auto const message = search::legacy_message(refused);
        return QString::fromLatin1(message.data(), static_cast<qsizetype>(message.size()));
    }
    return settings;
}

auto MainWindow::find_pair(FindPairRequest const& request) -> std::optional<QString> {
    if (not molecule_.has_value()) {
        return u"Please load a molecule"_s;
    }
    if (find_pair_running()) {
        return u"A Find Pair search is running"_s;
    }
    auto const settings = find_pair_settings(request);
    if (auto const* const error = std::get_if<QString>(&settings)) {
        return *error;
    }
    auto const& search_settings = std::get<search::FindPairSettings>(settings);
    show_find_pair(search_settings, search::legacy_find_pair(molecule_->sequence, search_settings,
                                                             defaults_.criteria, defaults_.weights));
    return std::nullopt;
}

// The search runs in another thread; the molecule cannot change meanwhile
// (load_molecule refuses), and the window waits for it when it closes
auto MainWindow::start_find_pair(FindPairRequest const& request) -> std::optional<QString> {
    if (not molecule_.has_value()) {
        return u"Please load a molecule"_s;
    }
    if (find_pair_running()) {
        return u"A Find Pair search is running"_s;
    }
    auto const settings = find_pair_settings(request);
    if (auto const* const error = std::get_if<QString>(&settings)) {
        return *error;
    }
    auto const& search_settings = std::get<search::FindPairSettings>(settings);
    pending_find_pair_ = search_settings;
    find_pair_interrupt_.store(false);
    find_pair_progress_.reset();
    find_pair_progress_.show();
    find_pair_watcher_.setFuture(QtConcurrent::run(
        [molecule = std::span<sequence::Nucleotide const>{molecule_->sequence}, search_settings,
         criteria = defaults_.criteria, weights = defaults_.weights, &interrupt = find_pair_interrupt_] {
            return search::legacy_find_pair(molecule, search_settings, criteria, weights, interrupt);
        }));
    return std::nullopt;
}

// until its result is shown: the thread may end before the window's
// event loop delivers the watcher's finished()
auto MainWindow::find_pair_running() const -> bool {
    return find_pair_watcher_.isRunning() or pending_find_pair_.has_value();
}

// PRIMER.exe's "Press key to interrupt search"
void MainWindow::cancel_find_pair() {
    find_pair_interrupt_.store(true);
}

void MainWindow::find_pair_finished() {
    find_pair_progress_.hide();
    if (not pending_find_pair_.has_value()) {
        return;
    }
    auto const settings = *pending_find_pair_;
    pending_find_pair_.reset();
    show_find_pair(settings, find_pair_watcher_.result());
}

void MainWindow::show_find_pair(search::FindPairSettings const& settings, search::FindPairResult const& result) {
    pairs_.show_pairs(result, defaults_.weights);
    report_.reset();  // PRIMER.exe prints nothing without a pair
    if (not result.pairs.empty() and molecule_.has_value()) {
        auto ranked = search::legacy_rank_pairs(result.pairs, defaults_.weights);
        auto const pair_count = ranked.size();
        auto best_pairs_text = [name = molecule_->name.toStdString(), size = molecule_->sequence.size(), settings,
                                ranked = std::move(ranked), criteria = defaults_.criteria](
                                   std::chrono::year_month_day const& date, std::size_t const best) {
            return report::legacy_find_pair_report({.molecule_name = name,
                                                    .molecule_size = size,
                                                    .settings = settings,
                                                    .best_pairs = report::legacy_printed_pairs(ranked, best)},
                                                   criteria, date);
        };
        // PRIMER.exe's "Print best: [N]": the best 20 at most
        report_ = ScreenReport{.file_name = u"FINDPAIR.PRN"_s,
                               .text = [best_pairs_text, pair_count](std::chrono::year_month_day const& date) {
                                   return best_pairs_text(date, report::legacy_print_best_default(pair_count));
                               },
                               .pair_count = pair_count,
                               .best_pairs_text = best_pairs_text};
    }
    pages_.setCurrentWidget(&pairs_);
}

void MainWindow::open_find_pair() {
    if (not molecule_.has_value() and not ask_for_molecule()) {
        return;
    }
    find_pair_.start(defaults_.find_pair_start_offset, defaults_.find_pair_end_offset);
}

void MainWindow::use_find_pair() {
    auto const request = find_pair_.request();
    if (not request.has_value()) {
        find_pair_.refuse(u"Invalid Position"_s);
        return;
    }
    if (auto const error = start_find_pair(*request); error.has_value()) {
        find_pair_.refuse(*error);
        return;
    }
    find_pair_.accept();
}

void MainWindow::open_find_enzymes() {
    if (not molecule_.has_value() and not ask_for_molecule()) {
        return;
    }
    find_enzymes_.start();
}

void MainWindow::use_find_enzymes() {
    if (auto const error = find_enzyme(find_enzymes_.enzyme(), find_enzymes_.mismatches()); error.has_value()) {
        find_enzymes_.refuse(*error);
        return;
    }
    find_enzymes_.accept();
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
    auto const complement = from_molecule_.strand() == sequence::Strand::complement;
    auto origin = u"%1 at %2%3"_s.arg(molecule_->name).arg(*position).arg(complement ? u" (complement)"_s : QString{});
    make_active(std::move(*primer), std::move(origin), sequence::PrimerLocation{.position = *position,
                                                             .length = defaults_.primer_length,
                                                             .strand = from_molecule_.strand()});
    if (active_primer_.has_value() and molecule_.has_value()) {  // always, here
        active_primer_->molecule_name = molecule_->name;
        active_primer_->molecule_path = QDir::toNativeSeparators(molecule_->path);
        active_primer_->position = static_cast<int>(*position);
        active_primer_->complement = complement;
    }
    start_picking();
}

// PRIMER.exe's MODIFY PRIMER POSITION screen for the active primer; while
// a second primer is picked, the one saved last to the Primer List is
// the partner of the A:B dimers line
void MainWindow::start_picking() {
    if (not active_primer_.has_value() or not active_primer_->location.has_value() or not molecule_.has_value()) {
        return;
    }
    auto partner = std::optional<PickingPartner>{};
    if (not primer_list_.empty()) {
        auto const& last = primer_list_.back();
        auto const bases = std::span{last.sequence}.subspan(last.extension_length);
        partner = PickingPartner{.primer = {bases.begin(), bases.end()},
                                 .name = u"primer %1 (%2)"_s.arg(primer_list_.size()).arg(letters(bases))};
    }
    auto const extension = std::span{active_primer_->sequence}.first(active_primer_->extension_length);
    picking_.start({.molecule_name = molecule_->name,
                    .molecule = molecule_->sequence,
                    .location = *active_primer_->location,
                    .criteria = active_primer_->criteria,
                    .partner = std::move(partner),
                    .extension = {extension.begin(), extension.end()}});
    pages_.setCurrentWidget(&picking_);
    picking_.setFocus();
}

auto MainWindow::modify_position() -> std::optional<QString> {
    if (active_primer_.has_value() and not active_primer_->original.empty()) {
        return u"Mutated primers cannot be moved"_s;  // PRIMER.exe's messages
    }
    if (active_primer_.has_value() and not active_primer_->location.has_value() and
        active_primer_->origin == u"Direct Entry"_s) {
        return u"Directly entered primers cannot be moved"_s;
    }
    if (not active_primer_.has_value() or not active_primer_->location.has_value() or not molecule_.has_value() or
        active_primer_->molecule_name != molecule_->name) {
        return u"Please create a primer from a molecule"_s;
    }
    start_picking();
    return std::nullopt;
}

void MainWindow::open_modify_position() {
    if (auto const error = modify_position(); error.has_value()) {
        QMessageBox::information(this, u"Modify Position"_s, *error);
    }
}

auto MainWindow::mutagenize() -> std::optional<QString> {
    if (not active_primer_.has_value() or not active_primer_->location.has_value() or not molecule_.has_value() or
        active_primer_->molecule_name != molecule_->name) {
        return u"Please create a primer from a molecule"_s;
    }
    auto const& primer = *active_primer_;
    if (not primer.original.empty()) {
        return u"Primer is already Mutated"_s;  // PRIMER.exe's message
    }
    // not recorded from PRIMER.exe yet (tests/legacy/README.md, "Mutagenize")
    if (primer.extension_length != 0) {
        return u"Mutagenize: primers without 5' extension only"_s;
    }
    mutagenize_source_ = primer;
    mutagenize_.start({.molecule_name = molecule_->name,
                       .molecule = molecule_->sequence,
                       .location = *primer.location,
                       .primer = primer.sequence,
                       .criteria = primer.criteria});
    pages_.setCurrentWidget(&mutagenize_);
    mutagenize_.setFocus();
    return std::nullopt;
}

auto MainWindow::view_molecule() -> std::optional<QString> {
    if (not molecule_.has_value()) {
        return u"Please load a molecule"_s;
    }
    molecule_view_.show_molecule({.name = molecule_->name,
                                  .sequence = molecule_->sequence,
                                  .file_text = molecule_->text,
                                  .topology = molecule_->topology});
    report_.reset();  // nothing to print on this screen (Print is not recorded)
    pages_.setCurrentWidget(&molecule_view_);
    return std::nullopt;
}

// as Find Enzymes, a molecule is asked for when none is loaded
void MainWindow::open_view_molecule() {
    if (not molecule_.has_value() and not ask_for_molecule()) {
        return;
    }
    static_cast<void>(view_molecule());  // a molecule is loaded: no error
}

auto MainWindow::remove_molecule() -> std::optional<QString> {
    if (not molecule_.has_value()) {
        return u"Molecule List is Empty"_s;  // PRIMER.exe's message
    }
    if (find_pair_running()) {
        return u"A Find Pair search is running"_s;  // it reads the molecule
    }
    molecule_.reset();
    molecule_label_.clear();
    // the screens that show the molecule close
    if (pages_.currentWidget() == &molecule_view_ or pages_.currentWidget() == &picking_ or
        pages_.currentWidget() == &mutagenize_) {
        analyze_active_primer();
    }
    return std::nullopt;
}

auto MainWindow::save_molecule_to(QString const& path) -> std::optional<QString> {
    if (not molecule_.has_value()) {
        return u"Molecule List is Empty"_s;  // PRIMER.exe's message
    }
    // a file of PRIMER.exe's own format is written back with its
    // features and description; another file as PRIMER.exe writes it
    // (tests/legacy/README.md, "Saving and creating molecules")
    auto native = io::parse_native_molecule(molecule_->text);
    if (not native.molecule.has_value()) {
        auto bases = std::string{};
        std::ranges::transform(molecule_->sequence, std::back_inserter(bases), sequence::to_letter);
        native.molecule = io::native_from_molecule(
            {.name = molecule_->name.toStdString(), .sequence = std::move(bases), .topology = molecule_->topology});
    }
    auto const text = io::write_native_molecule(*native.molecule);
    auto file = QFile{path};
    if (not file.open(QIODevice::WriteOnly) or
        file.write(text.data(), static_cast<qint64>(text.size())) != static_cast<qint64>(text.size())) {
        return u"cannot write the molecule to %1"_s.arg(path);
    }
    return std::nullopt;
}

// PRIMER.exe's SAVE/REMOVE MOLECULE screen: S saves the molecule (a file
// dialog, which asks before overwriting a file, as PRIMER.exe), R
// removes it
void MainWindow::open_save_remove() {
    if (not molecule_.has_value()) {
        QMessageBox::information(this, u"Save/Remove Molecule"_s, u"Molecule List is Empty"_s);
        return;
    }
    auto box = QMessageBox{QMessageBox::Question, u"Save/Remove Molecule"_s,
                           u"%1\n\nSave molecules to disk files and/or remove molecules from program list. "
                           u"Removing a molecule from the list unloads the molecule from program memory "
                           u"(it does not delete a disk file)."_s.arg(molecule_->name),
                           QMessageBox::Cancel, this};
    auto const* const save = box.addButton(u"&Save"_s, QMessageBox::AcceptRole);
    auto const* const remove = box.addButton(u"&Remove"_s, QMessageBox::DestructiveRole);
    box.exec();
    if (box.clickedButton() == save) {
        auto const path = QFileDialog::getSaveFileName(this, u"Save Molecule"_s, molecule_->path);
        if (path.isEmpty()) {
            return;
        }
        if (auto const error = save_molecule_to(path); error.has_value()) {
            QMessageBox::warning(this, u"Save/Remove Molecule"_s, *error);
        }
        return;
    }
    if (box.clickedButton() != remove) {
        return;
    }
    if (auto const error = remove_molecule(); error.has_value()) {
        QMessageBox::information(this, u"Save/Remove Molecule"_s, *error);
    }
}

void MainWindow::open_mutagenize() {
    if (not active_primer_.has_value()) {
        analyze_active_primer();  // asks for a primer
        return;
    }
    if (auto const error = mutagenize(); error.has_value()) {
        QMessageBox::information(this, u"Mutagenize"_s, *error);
    }
}

// S: as PRIMER.exe, the mutated primer is a new entry of the Primer
// List, and the active primer; the Mutagenize screen stays
void MainWindow::save_mutated_primer() {
    if (not mutagenize_source_.has_value()) {
        return;
    }
    active_primer_ = mutated_copy(*mutagenize_source_, mutagenize_.mutated());
    save_to_list();
    pages_.setCurrentWidget(&mutagenize_);
    mutagenize_.show_saved();
}

// I: the Primer Analysis of the mutated primer, whose report (Print to
// Disk) adds the mutagenesis information
void MainWindow::inspect_mutated_primer() {
    if (not mutagenize_source_.has_value()) {
        return;
    }
    active_primer_ = mutated_copy(*mutagenize_source_, mutagenize_.mutated());
    active_list_index_.reset();
    analyze_active_primer();
}

// the active primer follows the picking view's window; its 5'
// extension, if any, is kept
void MainWindow::use_picked_location() {
    if (not active_primer_.has_value() or not molecule_.has_value()) {
        return;
    }
    auto const location = picking_.location();
    auto const bases = sequence::primer_at(molecule_->sequence, location);
    if (not bases.has_value()) {
        return;
    }
    auto& primer = *active_primer_;
    active_list_index_.reset();  // a moved primer is a new one (S saves it as a new entry)
    primer.sequence.resize(primer.extension_length);
    primer.sequence.insert(primer.sequence.end(), bases->begin(), bases->end());
    primer.location = location;
    primer.position = static_cast<int>(location.position);
    primer.complement = location.strand == sequence::Strand::complement;
    primer.origin = u"%1 at %2%3"_s.arg(molecule_->name)
                        .arg(location.position)
                        .arg(primer.complement ? u" (complement)"_s : QString{});
}


// a new primer without extension, with the System Defaults' criteria
void MainWindow::make_active(std::vector<sequence::Nucleotide> sequence, QString origin,
                             std::optional<sequence::PrimerLocation> const location) {
    active_list_index_.reset();
    active_primer_ = ActivePrimer{.sequence = std::move(sequence),
                                  .extension_length = 0,
                                  .criteria = defaults_.criteria,
                                  .location = location,
                                  .origin = std::move(origin),
                                  .molecule_name = {},
                                  .molecule_path = {},
                                  .position = 1,
                                  .complement = false,
                                  .original = {}};
    analyze_active_primer();
}

// PRIMER.exe's S (Save+Exit): the active primer joins the Primer List
void MainWindow::save_to_list() {
    if (not active_primer_.has_value()) {
        analyze_active_primer();  // asks for a primer
        return;
    }
    primer_list_.push_back(*active_primer_);
    active_list_index_ = primer_list_.size() - 1;  // the active primer is that entry
    show_primer_list();
    if (pages_.currentWidget() == &picking_) {
        analyze_active_primer();  // as S on the picking screen
    }
}

void MainWindow::show_primer_list() {
    auto rows = QList<QList<QString>>{};
    for (auto const& primer : primer_list_) {
        auto const sequence = std::span{primer.sequence};
        auto const extension = letters(sequence.first(primer.extension_length));
        auto const bases = letters(sequence.subspan(primer.extension_length));
        rows.append({QString::number(rows.size() + 1), extension.isEmpty() ? bases : extension + u'-' + bases,
                     QString::number(primer.sequence.size()), primer.origin});
    }
    list_model_.set_rows(std::move(rows));
}

auto MainWindow::change_length(std::size_t const length) -> std::optional<QString> {
    if (active_primer_.has_value() and not active_primer_->original.empty()) {
        return u"Mutated primers cannot be moved"_s;  // not recorded for Length
    }
    if (not active_primer_.has_value() or not active_primer_->location.has_value() or not molecule_.has_value() or
        active_primer_->molecule_name != molecule_->name) {
        return u"Please create a primer from a molecule"_s;
    }
    auto& primer = *active_primer_;
    // the length includes the 5' extension, as the form says
    auto const bases = length > primer.extension_length ? length - primer.extension_length : 0;
    auto const [changed, message] = search::legacy_change_length(molecule_->sequence, *primer.location, bases);
    if (message == search::PickingMessage::primer_must_have_length) {
        return u"Primer must have length!"_s;  // PRIMER.exe's messages
    }
    if (message == search::PickingMessage::primer_too_long) {
        return u"Primer Too Long - Extends outside molecule"_s;
    }
    auto const primer_bases = sequence::primer_at(molecule_->sequence, changed);
    if (not primer_bases.has_value()) {
        return u"Primer Too Long - Extends outside molecule"_s;  // unreachable: the length was checked
    }
    primer.sequence.resize(primer.extension_length);
    primer.sequence.insert(primer.sequence.end(), primer_bases->begin(), primer_bases->end());
    primer.location = changed;
    // as PRIMER.exe, the Primer List entry of the active primer changes too
    if (active_list_index_.has_value() and *active_list_index_ < primer_list_.size()) {
        primer_list_.at(*active_list_index_) = primer;
        show_primer_list();
    }
    analyze_active_primer();
    return std::nullopt;
}

void MainWindow::open_change_length() {
    if (not active_primer_.has_value()) {
        analyze_active_primer();  // asks for a primer
        return;
    }
    auto accepted = false;
    constexpr auto minimum = 0;  // PRIMER.exe refuses it with a message
    constexpr auto maximum = 99;  // the form's two digits
    constexpr auto step = 1;
    auto const length = QInputDialog::getInt(this, u"Change Primer Length"_s, u"Primer length (incl 5' ext):"_s,
                                             static_cast<int>(active_primer_->sequence.size()), minimum, maximum,
                                             step, &accepted);
    if (not accepted) {
        return;
    }
    if (auto const error = change_length(static_cast<std::size_t>(length)); error.has_value()) {
        QMessageBox::information(this, u"Change Primer Length"_s, *error);
    }
}

void MainWindow::show_listed_primer(int const row) {
    if (row < 0 or std::cmp_greater_equal(row, primer_list_.size())) {
        return;
    }
    active_primer_ = primer_list_.at(static_cast<std::size_t>(row));
    active_list_index_ = static_cast<std::size_t>(row);
    analyze_active_primer();
}

auto MainWindow::analyze_pair(PairSelection const& selection) -> std::optional<QString> {
    if (selection.first >= primer_list_.size() or selection.second >= primer_list_.size()) {
        return u"Please save two primers to the Primer List"_s;
    }
    auto const& first = primer_list_.at(selection.first);
    auto const& second = primer_list_.at(selection.second);
    // the active primer's criteria, as the analysis screens
    auto const& criteria = active_primer_.has_value() ? active_primer_->criteria : defaults_.criteria;
    auto const primer_a = std::span{first.sequence}.subspan(first.extension_length);
    auto const primer_b = std::span{second.sequence}.subspan(second.extension_length);
    pair_.show_pair({.first = primer_a, .second = primer_b}, criteria);
    report_ = ScreenReport{.file_name = u"PAIR.PRN"_s,
                           .text = [primer_a = std::vector(primer_a.begin(), primer_a.end()),
                                    primer_b = std::vector(primer_b.begin(), primer_b.end()),
                                    criteria](std::chrono::year_month_day const& date) {
                               return report::legacy_pair_report({.first = primer_a, .second = primer_b}, criteria,
                                                                 date);
                           }};
    pages_.setCurrentWidget(&pair_);
    return std::nullopt;
}

void MainWindow::open_primer_pair() {
    constexpr auto pair_size = std::size_t{2};
    if (primer_list_.size() < pair_size) {
        QMessageBox::information(this, u"Primer Pair"_s, u"Please save two primers to the Primer List"_s);
        return;
    }
    pair_dialog_.start(list_texts());
}

// "N: SEQUENCE", a text per primer of the list
auto MainWindow::list_texts() const -> QList<QString> {
    auto texts = QList<QString>{};
    for (auto index = 0; index < list_model_.rowCount({}); ++index) {
        constexpr auto sequence_column = 1;
        texts.append(u"%1: %2"_s.arg(index + 1).arg(
            list_model_.data(list_model_.index(index, sequence_column), Qt::DisplayRole).toString()));
    }
    return texts;
}

auto MainWindow::analyze_multi(MultiSelection const& selection) -> std::optional<QString> {
    auto fields = analysis::MultiPrimerFields{};
    auto filled = std::size_t{0};
    for (auto index = std::size_t{0}; index < selection.size(); ++index) {
        auto const& choice = selection.at(index);
        if (not choice.has_value()) {
            continue;
        }
        if (*choice >= primer_list_.size()) {
            return u"Please save the primers to the Primer List"_s;
        }
        auto const& primer = primer_list_.at(*choice);
        fields.at(index) = std::span{primer.sequence}.subspan(primer.extension_length);
        ++filled;
    }
    constexpr auto min_filled = std::size_t{2};
    if (filled < min_filled) {
        return u"Please choose at least two primers"_s;
    }
    auto const& criteria = active_primer_.has_value() ? active_primer_->criteria : defaults_.criteria;
    multi_.show_primers(fields, criteria);
    // copies of the primers: the list may change before printing
    auto copies = std::array<std::optional<std::vector<sequence::Nucleotide>>, analysis::legacy_multi_primer_fields>{};
    std::ranges::transform(fields, copies.begin(), [](auto const& field) {
        return field.has_value() ? std::optional{std::vector(field->begin(), field->end())} : std::nullopt;
    });
    report_ = ScreenReport{.file_name = u"MULTI.PRN"_s,
                           .text = [copies = std::move(copies), criteria](std::chrono::year_month_day const& date) {
                               auto views = analysis::MultiPrimerFields{};
                               std::ranges::transform(copies, views.begin(), [](auto const& copy) {
                                   return copy.has_value()
                                              ? std::optional{std::span<sequence::Nucleotide const>{*copy}}
                                              : std::nullopt;
                               });
                               return report::legacy_multi_primer_report(views, criteria, date);
                           }};
    pages_.setCurrentWidget(&multi_);
    return std::nullopt;
}

auto MainWindow::show_glossary(QString const& name) -> bool {
    return glossary_.show_entry(name);
}

auto MainWindow::match_active_primer(MatchRequest const& request) -> std::optional<QString> {
    if (not active_primer_.has_value() or not active_primer_->location.has_value() or not molecule_.has_value()) {
        return u"Please create a primer from a molecule"_s;
    }
    auto const& location = *active_primer_->location;
    if (location.strand != sequence::Strand::top) {
        // primers B on the top strand: recorded (match_variations/), not
        // reproduced yet
        return u"Match Active Primer needs a primer A on the top strand"_s;
    }
    auto const primer_a = std::span{active_primer_->sequence}.subspan(active_primer_->extension_length);
    auto const settings = search::MatchSettings{.target = request.target,
                                                .start_offset = request.start_offset,
                                                .end_offset = request.end_offset,
                                                .strand = request.strand,
                                                .length = defaults_.primer_length,
                                                .gc_tolerance_percent = defaults_.gc_tolerance_percent,
                                                .tm_tolerance_celsius = defaults_.tm_tolerance_celsius};
    auto const list = search::legacy_match_active_primer({.molecule = molecule_->sequence, .primer_a = primer_a},
                                                         settings, active_primer_->criteria);
    match_.show_list(u"Primer A: 5'-%1-3', position %2 of %3; primers B (%7 strand) from %4 - %5 to %4 + %6"_s
                         .arg(letters(primer_a))
                         .arg(location.position)
                         .arg(molecule_->name)
                         .arg(request.target)
                         .arg(request.start_offset)
                         .arg(request.end_offset)
                         .arg(request.strand == sequence::Strand::top ? u"top"_s : u"complement"_s),
                     list);
    report_.reset();  // PRIMER.exe prints nothing without a primer B
    if (not list.empty()) {
        report_ = ScreenReport{
            .file_name = u"MATCH.PRN"_s,
            .text = [name = molecule_->name.toStdString(), size = molecule_->sequence.size(),
                     primer = std::vector(primer_a.begin(), primer_a.end()), position = location.position, settings,
                     criteria = active_primer_->criteria, list](std::chrono::year_month_day const& date) {
                return report::legacy_match_report({.molecule_name = name,
                                                    .molecule_size = size,
                                                    .primer_a = primer,
                                                    .primer_a_position = position,
                                                    .settings = settings,
                                                    .criteria = criteria,
                                                    .list = list},
                                                   date);
            }};
    }
    pages_.setCurrentWidget(&match_);
    return std::nullopt;
}

void MainWindow::open_match() {
    if (not active_primer_.has_value() or not active_primer_->location.has_value()) {
        QMessageBox::information(this, u"Match Active Primer"_s, u"Please create a primer from a molecule"_s);
        return;
    }
    match_dialog_.start(defaults_.match_start_offset, defaults_.match_end_offset);
}

void MainWindow::use_match() {
    auto const request = match_dialog_.request();
    if (not request.has_value()) {
        match_dialog_.refuse(u"Invalid Position"_s);
        return;
    }
    if (auto const error = match_active_primer(*request); error.has_value()) {
        match_dialog_.refuse(*error);
        return;
    }
    match_dialog_.accept();
}

void MainWindow::open_multi_primer() {
    constexpr auto min_primers = std::size_t{1};  // one primer may fill two fields
    if (primer_list_.size() < min_primers) {
        QMessageBox::information(this, u"Multi Primer"_s, u"Please save primers to the Primer List"_s);
        return;
    }
    multi_dialog_.start(list_texts());
}

void MainWindow::use_multi_primer() {
    if (auto const error = analyze_multi(multi_dialog_.selection()); error.has_value()) {
        QMessageBox::warning(this, u"Multi Primer"_s, *error);
    }
}

void MainWindow::use_primer_pair() {
    if (auto const error = analyze_pair(pair_dialog_.selection()); error.has_value()) {
        QMessageBox::warning(this, u"Primer Pair"_s, *error);
    }
}

auto MainWindow::homology(analysis::HomologyKind const kind) -> std::optional<QString> {
    if (not active_primer_.has_value() or not active_primer_->location.has_value() or not molecule_.has_value()) {
        return u"Please create a primer from a molecule"_s;
    }
    auto const primer = std::span{active_primer_->sequence}.subspan(active_primer_->extension_length);
    auto const settings = analysis::HomologySettings{
        .tm_filter = defaults_.exclude_tm_at_most_zero ? analysis::HomologyTmFilter::exclude_at_most_zero
                                                       : analysis::HomologyTmFilter::keep_all,
        .order = defaults_.sort_homologies_by_tm ? analysis::HomologyOrder::by_tm : analysis::HomologyOrder::by_position};
    auto const sites =
        analysis::legacy_homology_sites(primer, molecule_->sequence, *active_primer_->location, kind, settings);
    homology_.show_sites({.kind = kind, .primer = primer, .molecule_name = molecule_->name}, sites);
    // PRIMER.exe prints nothing without a site; the report of a
    // complement-strand primer is recorded, not reproduced yet
    report_.reset();
    if (not sites.empty() and active_primer_->location->strand == sequence::Strand::top) {
        report_ = ScreenReport{
            .file_name = kind == analysis::HomologyKind::uniqueness ? u"UNIQUE.PRN"_s : u"FALSEPRI.PRN"_s,
            .text = [name = molecule_->name.toStdString(), molecule = molecule_->sequence,
                     primer = std::vector(primer.begin(), primer.end()), position = active_primer_->location->position,
                     kind, tm_filter = settings.tm_filter, sites](std::chrono::year_month_day const& date) {
                return report::legacy_homology_report({.molecule_name = name,
                                                       .molecule_size = molecule.size(),
                                                       .primer = primer,
                                                       .position = position,
                                                       .kind = kind,
                                                       .tm_filter = tm_filter,
                                                       .sites = sites},
                                                      date);
            }};
    }
    pages_.setCurrentWidget(&homology_);
    return std::nullopt;
}

void MainWindow::show_homology(analysis::HomologyKind const kind) {
    if (auto const error = homology(kind); error.has_value()) {
        QMessageBox::information(this, u"Uniqueness and False Priming"_s, *error);
    }
}

void MainWindow::use_direct_entry() {
    make_active(direct_entry_.sequence(), u"Direct Entry"_s);
}

void MainWindow::modify_extension() {
    if (not active_primer_.has_value()) {
        analyze_active_primer();  // asks for a primer
        return;
    }
    auto const extension = std::span{active_primer_->sequence}.first(active_primer_->extension_length);
    auto text = QString{};
    for (auto const nucleotide : extension) {
        text.append(QChar::fromLatin1(sequence::to_letter(nucleotide)));
    }
    extension_.start(text);
}

// the primer keeps its criteria
void MainWindow::use_extension() {
    if (not active_primer_.has_value()) {
        return;
    }
    auto sequence = extension_.sequence();
    auto const extension_length = sequence.size();
    auto const primer = std::span{active_primer_->sequence}.subspan(active_primer_->extension_length);
    sequence.insert(sequence.end(), primer.begin(), primer.end());
    active_primer_->sequence = std::move(sequence);
    active_primer_->extension_length = extension_length;
    analyze_active_primer();
}

auto MainWindow::load_primer(QString const& path) -> std::optional<QString> {
    auto file = QFile{path};
    if (not file.open(QIODevice::ReadOnly)) {
        return u"cannot read %1"_s.arg(path);
    }
    auto const contents = file.readAll();
    auto bytes = std::vector<std::byte>{};
    bytes.reserve(static_cast<std::size_t>(contents.size()));
    std::ranges::transform(contents, std::back_inserter(bytes),
                           [](char const character) { return std::bit_cast<std::byte>(character); });
    auto const parsed = io::parse_pd2(bytes);
    auto const extension =
        parsed.primer.has_value() ? sequence::to_nucleotides(parsed.primer->extension) : std::nullopt;
    auto const primer = parsed.primer.has_value() ? sequence::to_nucleotides(parsed.primer->sequence) : std::nullopt;
    constexpr auto min_primer_length = std::size_t{2};  // the analyses need two bases
    if (not extension.has_value() or not primer.has_value() or primer->size() < min_primer_length) {
        return u"%1: not a primer file of Primer Designer 2.0 (.PD2)"_s.arg(path);
    }
    auto sequence = *extension;
    sequence.insert(sequence.end(), primer->begin(), primer->end());
    // the molecule of a .PD2 file may not be the loaded one: no location
    active_list_index_.reset();
    active_primer_ = ActivePrimer{.sequence = std::move(sequence),
                                  .extension_length = extension->size(),
                                  .criteria = parsed.primer->criteria,
                                  .location = std::nullopt,
                                  .origin = QFileInfo{path}.fileName(),
                                  .molecule_name = QString::fromStdString(parsed.primer->molecule_name),
                                  .molecule_path = QString::fromStdString(parsed.primer->molecule_path),
                                  .position = parsed.primer->position,
                                  .complement = parsed.primer->complement,
                                  .original = {}};
    analyze_active_primer();
    return std::nullopt;
}

void MainWindow::retrieve_primer() {
    auto const path = QFileDialog::getOpenFileName(this, u"Load Primer from Disk"_s, {},
                                                   u"Primer files (*.pd2 *.PD2);;All files (*)"_s);
    if (path.isEmpty()) {
        return;
    }
    if (auto const error = load_primer(path); error.has_value()) {
        QMessageBox::warning(this, u"Retrieve from disk"_s, *error);
    }
}

auto MainWindow::print_report_to(QString const& path) -> std::optional<QString> {
    if (not report_.has_value()) {
        return u"Nothing to print on this screen"_s;
    }
    return write_report(path, report_->text(today()));
}

auto MainWindow::print_report_to(QString const& path, std::size_t const best_pairs) -> std::optional<QString> {
    if (not report_.has_value()) {
        return u"Nothing to print on this screen"_s;
    }
    if (not report_->best_pairs_text or best_pairs == 0 or best_pairs > report::legacy_report_max_best_pairs) {
        return print_report_to(path);
    }
    return write_report(path, report_->best_pairs_text(today(), best_pairs));
}

auto MainWindow::save_primer_to(QString const& path) -> std::optional<QString> {
    if (not active_primer_.has_value()) {
        return u"Please create or retrieve a primer"_s;
    }
    auto const& active = *active_primer_;
    auto const sequence = std::span{active.sequence};
    auto const primer = io::Pd2Primer{.extension = letters(sequence.first(active.extension_length)).toStdString(),
                                .sequence = letters(sequence.subspan(active.extension_length)).toStdString(),
                                .molecule_name = active.molecule_name.toStdString(),
                                .molecule_path = active.molecule_path.toStdString(),
                                .position = active.position,
                                .complement = active.complement,
                                .criteria = active.criteria};
    auto const bytes = io::write_pd2(primer);
    auto contents = QByteArray{};
    contents.reserve(static_cast<qsizetype>(bytes.size()));
    for (auto const byte : bytes) {
        contents.append(std::bit_cast<char>(byte));
    }
    auto file = QFile{path};
    if (not file.open(QIODevice::WriteOnly) or file.write(contents) != contents.size()) {
        return u"cannot write the primer to %1"_s.arg(path);
    }
    return std::nullopt;
}

// PRIMER.exe asks for a file name and adds .PD2
void MainWindow::save_primer() {
    if (not active_primer_.has_value()) {
        analyze_active_primer();  // asks for a primer
        return;
    }
    auto path = QFileDialog::getSaveFileName(this, u"Save Active Primer"_s, {}, u"Primer files (*.PD2)"_s);
    if (path.isEmpty()) {
        return;
    }
    if (QFileInfo{path}.suffix().isEmpty()) {
        path += u".PD2"_s;
    }
    if (auto const error = save_primer_to(path); error.has_value()) {
        QMessageBox::warning(this, u"Save Active Primer"_s, *error);
    }
}

// Print Control -> Print to Disk: PRIMER.exe writes its report to a
// text file
void MainWindow::print_report() {
    if (not report_.has_value()) {
        QMessageBox::information(this, u"Print to Disk"_s, u"Nothing to print on this screen"_s);
        return;
    }
    // the Find Pair list: PRIMER.exe's "N Pairs found.  Print best: [M]"
    auto best_pairs = std::size_t{0};
    if (report_->best_pairs_text) {
        auto accepted = false;
        constexpr auto minimum = 1;
        constexpr auto step = 1;
        auto const answer = QInputDialog::getInt(
            this, u"Print to Disk"_s, u"%1 Pairs found.  Print best:"_s.arg(report_->pair_count),
            static_cast<int>(report::legacy_print_best_default(report_->pair_count)), minimum,
            static_cast<int>(report::legacy_report_max_best_pairs), step, &accepted);
        if (not accepted) {
            return;
        }
        best_pairs = static_cast<std::size_t>(answer);
    }
    auto const path = QFileDialog::getSaveFileName(this, u"Print to Disk"_s, report_->file_name,
                                                   u"Reports (*.prn *.PRN *.txt);;All files (*)"_s);
    if (path.isEmpty()) {
        return;
    }
    if (auto const error = print_report_to(path, best_pairs); error.has_value()) {
        QMessageBox::warning(this, u"Print to Disk"_s, *error);
    }
}

void MainWindow::analyze_active_primer() {
    if (not active_primer_.has_value()) {
        report_.reset();
        pages_.setCurrentWidget(&no_primer_);
        return;
    }
    analysis_.show_primer({.sequence = active_primer_->sequence, .extension_length = active_primer_->extension_length},
                          active_primer_->criteria);
    if (not active_primer_->original.empty() and molecule_.has_value() and active_primer_->location.has_value()) {
        report_ = ScreenReport{
            .file_name = u"INSPECT.PRN"_s,
            .text = [primer = *active_primer_, name = molecule_->name.toStdString(),
                     size = molecule_->sequence.size()](std::chrono::year_month_day const& date) {
                return report::legacy_mutated_inspect_report({.original = primer.original,
                                                              .mutated = primer.sequence,
                                                              .molecule_name = name,
                                                              .molecule_size = size,
                                                              .position = primer.location->position},
                                                             primer.criteria, date);
            }};
        pages_.setCurrentWidget(&analysis_);
        return;
    }
    report_ = ScreenReport{.file_name = u"INSPECT.PRN"_s,
                           .text = [primer = *active_primer_](std::chrono::year_month_day const& date) {
                               return report::legacy_inspect_report(
                                   {.sequence = primer.sequence, .extension_length = primer.extension_length},
                                   primer.criteria, date);
                           }};
    pages_.setCurrentWidget(&analysis_);
}

}  // namespace primerpicker::gui
