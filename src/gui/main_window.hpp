#pragma once

#include "core/analysis/legacy_verdict.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/topology.hpp"
#include "gui/analysis_view.hpp"
#include "gui/colour_scheme.hpp"
#include "gui/from_molecule_dialog.hpp"
#include "gui/sequence_dialog.hpp"

#include <QActionGroup>
#include <QLabel>
#include <QMainWindow>
#include <QMenu>
#include <QObject>
#include <QSettings>
#include <QStackedWidget>
#include <QString>

#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace primerpicker::gui {

// The primer the analyses work on: the whole sequence 5' to 3', whose
// first extension_length nucleotides are the 5' extension, and its
// criteria. As in a .PD2 file, a primer has its own criteria (Modify ->
// Active Primer: Criteria): the System Defaults' when it is made, the
// saved ones when it is retrieved from disk.
struct ActivePrimer {
    std::vector<sequence::Nucleotide> sequence;
    std::size_t extension_length = 0;
    analysis::LegacyCriteria criteria;
};

// A molecule (template) loaded from a file
struct LoadedMolecule {
    QString name;
    std::vector<sequence::Nucleotide> sequence;
    sequence::Topology topology = sequence::Topology::linear;
};

// The main window: Primer Designer's main menu as a menu bar. Items
// that are not implemented yet are disabled. With settings, the window
// keeps its own preferences there (the colour scheme); Primer
// Designer's settings come from a pd2.cfg file (load_settings).
class MainWindow final : public QMainWindow {
    Q_OBJECT

public:
    // without settings: preferences are not kept
    MainWindow();
    explicit MainWindow(QSettings& settings);

    // Loads a molecule file (GenBank, FASTA or plain sequence) as the
    // molecule primers are made from; the error, if it was not loaded
    [[nodiscard]] auto load_molecule(QString const& path) -> std::optional<QString>;

    // Loads Primer Designer's settings from a pd2.cfg file (criteria,
    // primer length, tolerances, search ranges and weights); the error,
    // if it was not loaded (the settings are then kept)
    [[nodiscard]] auto load_settings(QString const& path) -> std::optional<QString>;

    // Retrieves a primer saved by Primer Designer (.PD2 file) as the
    // active primer, and shows its analysis; the error, if it was not
    // loaded (the active primer is then kept)
    [[nodiscard]] auto load_primer(QString const& path) -> std::optional<QString>;

    // Writes PRIMER.exe's One Primer report of the active primer (Print
    // to Disk), dated today; the error, if it was not written
    [[nodiscard]] auto print_report_to(QString const& path) -> std::optional<QString>;

private:
    void add_legacy_menus();
    void add_colour_scheme_menu(QMenu& parent);
    void choose_colour_scheme(ColourScheme scheme);
    void use_direct_entry();
    void modify_extension();
    void use_extension();
    void retrieve_primer();
    void print_report();
    void make_active(std::vector<sequence::Nucleotide> sequence, std::size_t extension_length);
    auto ask_for_molecule() -> bool;
    void create_from_molecule();
    void use_from_molecule();
    void analyze_active_primer();

    std::optional<std::reference_wrapper<QSettings>> settings_;
    io::LegacySystemDefaults defaults_;
    std::optional<ActivePrimer> active_primer_;
    std::optional<LoadedMolecule> molecule_;
    // members are destroyed in reverse order: the pages before the
    // stack that holds them
    QStackedWidget pages_;
    QLabel no_primer_;
    AnalysisView analysis_;
    QLabel molecule_label_;
    QLabel settings_label_;
    QActionGroup colour_schemes_;
    SequenceDialog direct_entry_;
    SequenceDialog extension_;
    FromMoleculeDialog from_molecule_;
};

}  // namespace primerpicker::gui
