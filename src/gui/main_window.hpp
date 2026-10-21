#pragma once

#include "core/analysis/legacy_homology.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/enzymes/site_search.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "core/sequence/topology.hpp"
#include "gui/analysis_view.hpp"
#include "gui/colour_scheme.hpp"
#include "gui/criteria_dialog.hpp"
#include "gui/find_enzymes_dialog.hpp"
#include "gui/find_pair_dialog.hpp"
#include "gui/find_pair_view.hpp"
#include "gui/from_molecule_dialog.hpp"
#include "gui/glossary_dialog.hpp"
#include "gui/homology_view.hpp"
#include "gui/match_dialog.hpp"
#include "gui/match_view.hpp"
#include "gui/multi_primer_dialog.hpp"
#include "gui/multi_primer_view.hpp"
#include "gui/pair_dialog.hpp"
#include "gui/pair_view.hpp"
#include "gui/sequence_dialog.hpp"
#include "gui/sites_view.hpp"
#include "gui/system_defaults_dialog.hpp"
#include "gui/text_table_model.hpp"

#include <QActionGroup>
#include <QDockWidget>
#include <QLabel>
#include <QList>
#include <QMainWindow>
#include <QMenu>
#include <QObject>
#include <QSettings>
#include <QStackedWidget>
#include <QTableView>
#include <QString>

#include <cstddef>
#include <chrono>
#include <functional>
#include <string>
#include <optional>
#include <vector>

namespace primerpicker::gui {

// The primer the analyses work on: the whole sequence 5' to 3', whose
// first extension_length nucleotides are the 5' extension, and its
// criteria. As in a .PD2 file, a primer has its own criteria (Modify ->
// Active Primer: Criteria): the System Defaults' when it is made, the
// saved ones when it is retrieved from disk. A primer made from the
// loaded molecule has its location there (without its extension).
struct ActivePrimer {
    std::vector<sequence::Nucleotide> sequence;
    std::size_t extension_length = 0;
    analysis::LegacyCriteria criteria;
    std::optional<sequence::PrimerLocation> location;
    QString origin;  // shown in the Primer List: Direct Entry, a molecule and position, a file
    // what a .PD2 file keeps: the source molecule (none for Direct
    // Entry), and the primer's position (1 for Direct Entry) and strand
    QString molecule_name;
    QString molecule_path;
    int position = 1;
    bool complement = false;
};

// A molecule (template) loaded from a file
struct LoadedMolecule {
    QString name;
    QString path;  // its file
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

    // Writes the System Defaults as a pd2.cfg file of Primer Designer;
    // the error, if it was not written
    [[nodiscard]] auto save_settings_to(QString const& path) const -> std::optional<QString>;

    // Retrieves a primer saved by Primer Designer (.PD2 file) as the
    // active primer, and shows its analysis; the error, if it was not
    // loaded (the active primer is then kept)
    [[nodiscard]] auto load_primer(QString const& path) -> std::optional<QString>;

    // Print to Disk: writes PRIMER.exe's report of the current screen
    // (One Primer, Primer Pair, Multiple Primers, Find Pair, Uniqueness,
    // False Priming, Match Active Primer), dated today; the error, if it
    // was not written (a screen without report, or nothing to print, as
    // PRIMER.exe for an empty list)
    [[nodiscard]] auto print_report_to(QString const& path) -> std::optional<QString>;

    // Primer/Create -> Save Active Primer: writes the active primer as a
    // .PD2 file of Primer Designer; the error, if it was not written
    [[nodiscard]] auto save_primer_to(QString const& path) -> std::optional<QString>;

    // Molecule-Ops -> Find Enzymes on the loaded molecule: shows the true
    // or almost sites of an enzyme; the error, if there was no search
    [[nodiscard]] auto find_enzyme(QString const& name, enzymes::Mismatches mismatches) -> std::optional<QString>;

    // Primer/Create -> Find Pair on the loaded molecule, with the
    // settings' criteria, primer length, tolerances and weights: shows
    // the ranked pairs; the error, if there was no search
    [[nodiscard]] auto find_pair(FindPairRequest const& request) -> std::optional<QString>;

    // The Uniqueness or False Priming screen of the active primer, made
    // from the loaded molecule (its 5' extension aside), with the
    // settings' homology options; the error, if there was no search
    [[nodiscard]] auto homology(analysis::HomologyKind kind) -> std::optional<QString>;

    // Analyze -> Primer Pair of two primers of the Primer List (their 5'
    // extensions aside); the error, if there was no analysis
    [[nodiscard]] auto analyze_pair(PairSelection const& selection) -> std::optional<QString>;

    // Analyze -> Multi Primer of two to four primers of the Primer List
    // (their 5' extensions aside); the error, if there was no analysis
    [[nodiscard]] auto analyze_multi(MultiSelection const& selection) -> std::optional<QString>;

    // Primer/Create -> Match Active Primer: the primers B (complement
    // strand) that match the active primer, made from the loaded
    // molecule on its top strand (5' extension aside), with the
    // settings' criteria, length and tolerances; the error, if there was
    // no search
    [[nodiscard]] auto match_active_primer(MatchRequest const& request) -> std::optional<QString>;

    // Help -> Glossary at the built-in definition NAME (a key or an
    // alias; the first entry if there is none); whether NAME was found
    auto show_glossary(QString const& name) -> bool;

private:
    void add_legacy_menus();
    void add_colour_scheme_menu(QMenu& parent);
    void choose_colour_scheme(ColourScheme scheme);
    void use_direct_entry();
    void modify_extension();
    void use_extension();
    void retrieve_primer();
    void save_primer();
    void modify_criteria();
    void use_criteria();
    void open_system_defaults();
    void use_system_defaults();
    void load_settings_file();
    void save_settings_file();
    void open_find_enzymes();
    void use_find_enzymes();
    void open_find_pair();
    void use_find_pair();
    void print_report();
    void make_active(std::vector<sequence::Nucleotide> sequence, QString origin,
                     std::optional<sequence::PrimerLocation> location = std::nullopt);
    void save_to_list();
    void show_listed_primer(int row);
    void open_primer_pair();
    void use_primer_pair();
    void open_multi_primer();
    void use_multi_primer();
    void open_match();
    void use_match();
    [[nodiscard]] auto list_texts() const -> QList<QString>;
    void show_homology(analysis::HomologyKind kind);
    auto ask_for_molecule() -> bool;
    void create_from_molecule();
    void use_from_molecule();
    void analyze_active_primer();

    std::optional<std::reference_wrapper<QSettings>> settings_;
    io::LegacySystemDefaults defaults_;
    QString settings_name_;  // the settings file's name; empty: built-in
    std::optional<ActivePrimer> active_primer_;
    std::optional<LoadedMolecule> molecule_;
    // the report of the current screen: its usual file name, and its text
    // for a date
    struct ScreenReport {
        QString file_name;
        std::function<std::string(std::chrono::year_month_day const&)> text;
    };
    std::optional<ScreenReport> report_;
    std::vector<ActivePrimer> primer_list_;  // the Primer List, in saving order
    // members are destroyed in reverse order: the pages before the
    // stack that holds them
    QStackedWidget pages_;
    QLabel no_primer_;
    AnalysisView analysis_;
    SitesView sites_;
    FindPairView pairs_;
    HomologyView homology_;
    PairView pair_;
    MultiPrimerView multi_;
    MatchView match_;
    // the Primer List panel: the dock before the view it holds
    QDockWidget list_dock_;
    TextTableModel list_model_;
    QTableView list_view_;
    QLabel molecule_label_;
    QLabel settings_label_;
    QActionGroup colour_schemes_;
    SequenceDialog direct_entry_;
    SequenceDialog extension_;
    FromMoleculeDialog from_molecule_;
    FindEnzymesDialog find_enzymes_;
    FindPairDialog find_pair_;
    PairDialog pair_dialog_;
    MultiPrimerDialog multi_dialog_;
    MatchDialog match_dialog_;
    GlossaryDialog glossary_;
    CriteriaDialog criteria_dialog_;
    SystemDefaultsDialog defaults_dialog_;
};

}  // namespace primerpicker::gui
