#include "gui/main_window.hpp"

#include "core/analysis/legacy_homology.hpp"
#include "core/enzymes/rebase_enzymes.hpp"
#include "core/enzymes/site_search.hpp"
#include "core/io/molecule_file.hpp"
#include "core/search/legacy_find_pair.hpp"
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
#include "gui/sequence_dialog.hpp"

#include <QAction>
#include <QAnyStringView>
#include <QChar>
#include <QDate>
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
#include <chrono>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <optional>
#include <span>
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

MainWindow::MainWindow()
    : colour_schemes_{this},
      direct_entry_{*this, direct_entry_kind},
      extension_{*this, extension_kind},
      from_molecule_{*this},
      find_enzymes_{*this},
      find_pair_{*this} {
    setWindowTitle(u"PrimerPicker"_s);
    add_legacy_menus();
    no_primer_.setObjectName(u"no_primer_label"_s);
    no_primer_.setText(QString::fromUtf16(no_primer_message));
    no_primer_.setAlignment(Qt::AlignCenter);
    pages_.addWidget(&no_primer_);
    pages_.addWidget(&analysis_);
    pages_.addWidget(&sites_);
    pages_.addWidget(&pairs_);
    pages_.addWidget(&homology_);
    setCentralWidget(&pages_);
    molecule_label_.setObjectName(u"molecule_label"_s);
    statusBar()->addWidget(&molecule_label_);
    settings_label_.setObjectName(u"settings_label"_s);
    statusBar()->addPermanentWidget(&settings_label_);
    QObject::connect(&direct_entry_, &QDialog::accepted, this, &MainWindow::use_direct_entry);
    QObject::connect(&extension_, &QDialog::accepted, this, &MainWindow::use_extension);
    QObject::connect(&find_enzymes_, &FindEnzymesDialog::enzyme_entered, this, &MainWindow::use_find_enzymes);
    QObject::connect(&find_pair_, &FindPairDialog::search_entered, this, &MainWindow::use_find_pair);
    QObject::connect(&analysis_, &AnalysisView::print_requested, this, &MainWindow::print_report);
    QObject::connect(&analysis_, &AnalysisView::uniqueness_requested, this,
                     [this] { show_homology(analysis::HomologyKind::uniqueness); });
    QObject::connect(&analysis_, &AnalysisView::false_priming_requested, this,
                     [this] { show_homology(analysis::HomologyKind::false_priming); });
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
    QObject::connect(&direct_entry, &QAction::triggered, &direct_entry_, [this] { direct_entry_.start(); });
    add_disabled_items(primer_create, {u"Match &Active Primer"_s});
    auto const& find_pair_item = add_item(primer_create, {.text = u"Find &Pair...", .name = u"find_pair"});
    QObject::connect(&find_pair_item, &QAction::triggered, this, &MainWindow::open_find_pair);
    add_disabled_items(primer_create, {u"&Save Active Primer"_s});
    auto const& retrieve = add_item(primer_create, {.text = u"&Retrieve from disk...", .name = u"retrieve_primer"});
    QObject::connect(&retrieve, &QAction::triggered, this, &MainWindow::retrieve_primer);

    auto& modify = *menu_bar.addMenu(u"&Modify"_s);
    modify.addSection(u"Active Primer:"_s);
    add_disabled_items(modify, {u"&Position"_s, u"&Length"_s});
    auto const& extension = add_item(modify, {.text = u"5'-&extension...", .name = u"modify_extension"});
    QObject::connect(&extension, &QAction::triggered, this, &MainWindow::modify_extension);
    add_disabled_items(modify, {u"&Criteria"_s, u"M&utagenize"_s});
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
    auto const& find_enzymes = add_item(molecule_ops, {.text = u"Find &Enzymes...", .name = u"find_enzymes"});
    QObject::connect(&find_enzymes, &QAction::triggered, this, &MainWindow::open_find_enzymes);
    add_disabled_items(molecule_ops, {u"&View/Edit"_s, u"&Save/Remove"_s});

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
    pages_.setCurrentWidget(&sites_);
    return std::nullopt;
}

auto MainWindow::find_pair(FindPairRequest const& request) -> std::optional<QString> {
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
    auto const result = search::legacy_find_pair(molecule_->sequence, settings, defaults_.criteria);
    pairs_.show_pairs(result, defaults_.weights);
    pages_.setCurrentWidget(&pairs_);
    return std::nullopt;
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
    if (auto const error = find_pair(*request); error.has_value()) {
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
    make_active(std::move(*primer), sequence::PrimerLocation{.position = *position,
                                                             .length = defaults_.primer_length,
                                                             .strand = from_molecule_.strand()});
}

// a new primer without extension, with the System Defaults' criteria
void MainWindow::make_active(std::vector<sequence::Nucleotide> sequence,
                             std::optional<sequence::PrimerLocation> const location) {
    active_primer_ = ActivePrimer{
        .sequence = std::move(sequence), .extension_length = 0, .criteria = defaults_.criteria, .location = location};
    analyze_active_primer();
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
    pages_.setCurrentWidget(&homology_);
    return std::nullopt;
}

void MainWindow::show_homology(analysis::HomologyKind const kind) {
    if (auto const error = homology(kind); error.has_value()) {
        QMessageBox::information(this, u"Uniqueness and False Priming"_s, *error);
    }
}

void MainWindow::use_direct_entry() {
    make_active(direct_entry_.sequence());
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
    active_primer_ = ActivePrimer{.sequence = std::move(sequence),
                                  .extension_length = extension->size(),
                                  .criteria = parsed.primer->criteria,
                                  .location = std::nullopt};
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
    if (not active_primer_.has_value()) {
        return u"Please create or retrieve a primer"_s;
    }
    auto const today = QDate::currentDate();
    auto const date = std::chrono::year_month_day{std::chrono::year{today.year()},
                                                  std::chrono::month{static_cast<unsigned>(today.month())},
                                                  std::chrono::day{static_cast<unsigned>(today.day())}};
    auto const report = report::legacy_inspect_report(
        {.sequence = active_primer_->sequence, .extension_length = active_primer_->extension_length},
        active_primer_->criteria, date);
    auto file = QFile{path};
    if (not file.open(QIODevice::WriteOnly) or
        file.write(report.data(), static_cast<qint64>(report.size())) != static_cast<qint64>(report.size())) {
        return u"cannot write the report to %1"_s.arg(path);
    }
    return std::nullopt;
}

// Print Control -> Print to Disk: PRIMER.exe writes its report to a
// text file
void MainWindow::print_report() {
    auto const path = QFileDialog::getSaveFileName(this, u"Print to Disk"_s, u"INSPECT.PRN"_s,
                                                   u"Reports (*.prn *.PRN *.txt);;All files (*)"_s);
    if (path.isEmpty()) {
        return;
    }
    if (auto const error = print_report_to(path); error.has_value()) {
        QMessageBox::warning(this, u"Print to Disk"_s, *error);
    }
}

void MainWindow::analyze_active_primer() {
    if (not active_primer_.has_value()) {
        pages_.setCurrentWidget(&no_primer_);
        return;
    }
    analysis_.show_primer({.sequence = active_primer_->sequence, .extension_length = active_primer_->extension_length},
                          active_primer_->criteria);
    pages_.setCurrentWidget(&analysis_);
}

}  // namespace primerpicker::gui
