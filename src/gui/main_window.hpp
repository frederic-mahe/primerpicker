#pragma once

#include "core/analysis/legacy_homology.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/enzymes/site_search.hpp"
#include "core/io/pd2_cfg.hpp"
#include "core/search/legacy_find_pair.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "core/sequence/topology.hpp"
#include "gui/about_dialog.hpp"
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
#include "gui/molecule_view.hpp"
#include "gui/multi_primer_view.hpp"
#include "gui/mutagenize_view.hpp"
#include "gui/pair_dialog.hpp"
#include "gui/pair_view.hpp"
#include "gui/picking_view.hpp"
#include "gui/return_key_filter.hpp"
#include "gui/sequence_dialog.hpp"
#include "gui/sites_view.hpp"
#include "gui/system_defaults_dialog.hpp"
#include "gui/text_table_model.hpp"

#include <QActionGroup>
#include <QByteArray>
#include <QCloseEvent>
#include <QDockWidget>
#include <QFutureWatcher>
#include <QLabel>
#include <QList>
#include <QMainWindow>
#include <QMenu>
#include <QProgressDialog>
#include <QObject>
#include <QSettings>
#include <QStackedWidget>
#include <QTableView>
#include <QString>

#include <atomic>
#include <cstddef>
#include <chrono>
#include <functional>
#include <string>
#include <optional>
#include <variant>
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
    // a mutated primer (saved on the Mutagenize screen): the primer it
    // was made from; empty for other primers
    std::vector<sequence::Nucleotide> original;
};

// A molecule (template) loaded from a file
struct LoadedMolecule {
    QString name;
    QString path;  // its file
    std::string text;  // the file's contents (View/Edit: features, annotations)
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
    ~MainWindow() override;

    MainWindow(MainWindow const&) = delete;
    MainWindow(MainWindow&&) = delete;
    auto operator=(MainWindow const&) -> MainWindow& = delete;
    auto operator=(MainWindow&&) -> MainWindow& = delete;

    // Loads a molecule file (GenBank, FASTA or plain sequence) as the
    // molecule primers are made from; the error, if it was not loaded
    [[nodiscard]] auto load_molecule(QString const& path) -> std::optional<QString>;

    // Help -> Example Molecule (not in PRIMER.exe): loads pUC19, built
    // into the program (data/molecules/pUC19.gb), as if from a file
    // named pUC19.gb; the error, if it was not loaded
    [[nodiscard]] auto load_example_molecule() -> std::optional<QString>;

    // Loads Primer Designer's settings from a pd2.cfg file (criteria,
    // primer length, tolerances, search ranges and weights); the error,
    // if it was not loaded (the settings are then kept)
    [[nodiscard]] auto load_settings(QString const& path) -> std::optional<QString>;

    // Writes the System Defaults as a pd2.cfg file of Primer Designer;
    // the error, if it was not written
    [[nodiscard]] auto save_settings_to(QString const& path) const -> std::optional<QString>;

    // Whether the System Defaults differ from those loaded or saved last
    // (the built-in ones at first): quitting then asks whether to save
    // them, as PRIMER.exe
    [[nodiscard]] auto settings_changed() const -> bool;

    // Where quitting offers to save changed System Defaults: PD2.CFG, in
    // the folder of the settings file loaded last, if any
    [[nodiscard]] auto settings_save_path() const -> QString;

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
    // The same; on the Find Pair list, with the answer to PRIMER.exe's
    // "N Pairs found.  Print best: [M]" (1 to 9999) instead of M
    [[nodiscard]] auto print_report_to(QString const& path, std::size_t best_pairs) -> std::optional<QString>;

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

    // The same search in another thread, behind a progress dialog whose
    // Cancel interrupts it (the form's OK does this); the error, if it
    // did not start. The result is shown when the search ends.
    [[nodiscard]] auto start_find_pair(FindPairRequest const& request) -> std::optional<QString>;
    [[nodiscard]] auto find_pair_running() const -> bool;
    void cancel_find_pair();

    // The Uniqueness or False Priming screen of the active primer, made
    // from the loaded molecule (its 5' extension aside), with the
    // settings' homology options; the error, if there was no search
    [[nodiscard]] auto homology(analysis::HomologyKind kind) -> std::optional<QString>;

    // Analyze -> Primer Pair of two primers of the Primer List (their 5'
    // extensions included, as PRIMER.exe); the error, if there was no
    // analysis
    [[nodiscard]] auto analyze_pair(PairSelection const& selection) -> std::optional<QString>;

    // Analyze -> Multi Primer of two to four primers of the Primer List
    // (their 5' extensions included, as PRIMER.exe); the error, if there
    // was no analysis
    [[nodiscard]] auto analyze_multi(MultiSelection const& selection) -> std::optional<QString>;

    // Primer/Create -> Match Active Primer: the primers B (complement
    // strand) that match the active primer, made from the loaded
    // molecule on its top strand (5' extension aside), with the
    // settings' criteria, length and tolerances; the error, if there was
    // no search
    [[nodiscard]] auto match_active_primer(MatchRequest const& request) -> std::optional<QString>;

    // Modify -> Position: the picking view for the active primer, made
    // from the loaded molecule (From Molecule opens it too); the error,
    // if there was none
    [[nodiscard]] auto modify_position() -> std::optional<QString>;

    // Modify -> Length: the active primer's length (5' extension
    // included), its 5' end kept; the Primer List entry it comes from
    // changes too, as in PRIMER.exe; the error, if it was not changed
    [[nodiscard]] auto change_length(std::size_t length) -> std::optional<QString>;

    // Modify -> Mutagenize: the Mutagenize Primer screen for the active
    // primer, made from the loaded molecule on its top strand, without
    // 5' extension; the error, if there was none
    [[nodiscard]] auto mutagenize() -> std::optional<QString>;

    // Molecule-Ops -> View/Edit: the loaded molecule's views (read only);
    // the error, if there is no molecule
    [[nodiscard]] auto view_molecule() -> std::optional<QString>;

    // Molecule-Ops -> Save/Remove, R: unloads the loaded molecule (its
    // file stays); the error, if none is loaded ("Molecule List is
    // Empty", as PRIMER.exe)
    [[nodiscard]] auto remove_molecule() -> std::optional<QString>;

    // Molecule-Ops -> Save/Remove, S: writes the loaded molecule to PATH
    // in PRIMER.exe's own format, as PRIMER.exe does; the error, if none
    // is loaded ("Molecule List is Empty") or if the file cannot be
    // written
    [[nodiscard]] auto save_molecule_to(QString const& path) -> std::optional<QString>;

    // Help -> Glossary at the built-in definition NAME (a key or an
    // alias; the first entry if there is none); whether NAME was found
    auto show_glossary(QString const& name) -> bool;

protected:
    // PRIMER.exe's "System Defaults Changed But Not Saved! SAVE CURRENT
    // DEFAULTS (Y/N) ?": Yes opens a save dialog, No quits, Cancel stays
    void closeEvent(QCloseEvent* event) override;

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
    void find_pair_finished();
    void show_find_pair(search::FindPairSettings const& settings, search::FindPairResult const& result);
    [[nodiscard]] auto find_pair_settings(FindPairRequest const& request) const
        -> std::variant<search::FindPairSettings, QString>;
    void print_report();
    void make_active(std::vector<sequence::Nucleotide> sequence, QString origin,
                     std::optional<sequence::PrimerLocation> location = std::nullopt);
    void save_to_list();
    void show_primer_list();
    void open_change_length();
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
    void start_picking();
    void open_modify_position();
    void use_picked_location();
    void open_mutagenize();
    void open_view_molecule();
    // A molecule file's contents, read from DISPLAY_PATH (in messages),
    // kept with PATH (the default of Save, primer files)
    struct MoleculeText {
        QByteArray bytes;
        QString display_path;
        QString path;
    };
    auto use_molecule_text(MoleculeText const& text) -> std::optional<QString>;
    void open_save_remove();
    void save_mutated_primer();
    void inspect_mutated_primer();
    void analyze_active_primer();

    std::optional<std::reference_wrapper<QSettings>> settings_;
    io::LegacySystemDefaults defaults_;
    QString settings_name_;  // the settings file's name; empty: built-in
    QString settings_path_;  // the settings file loaded last; empty: none
    io::LegacySystemDefaults saved_defaults_;  // the System Defaults loaded or saved last
    std::optional<ActivePrimer> active_primer_;
    std::optional<ActivePrimer> mutagenize_source_;  // the primer the Mutagenize screen was opened on
    std::optional<LoadedMolecule> molecule_;
    // the report of the current screen: its usual file name, and its text
    // for a date
    struct ScreenReport {
        QString file_name;
        std::function<std::string(std::chrono::year_month_day const&)> text;
        // the Find Pair list: the pairs found, and the report of the best N
        std::size_t pair_count = 0;
        std::function<std::string(std::chrono::year_month_day const&, std::size_t)> best_pairs_text = {};
    };
    std::optional<ScreenReport> report_;
    std::vector<ActivePrimer> primer_list_;  // the Primer List, in saving order
    std::optional<std::size_t> active_list_index_;  // the entry the active primer is, if any
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
    PickingView picking_;
    MutagenizeView mutagenize_;
    MoleculeView molecule_view_;
    // the Primer List panel: the dock before the view it holds
    QDockWidget list_dock_;
    TextTableModel list_model_;
    QTableView list_view_;
    // Return and Enter in the list pick the primer (as activated() does
    // on Windows and Linux, but not on macOS)
    ReturnKeyFilter list_return_;
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
    // the Find Pair search in another thread
    std::atomic<bool> find_pair_interrupt_{false};
    std::optional<search::FindPairSettings> pending_find_pair_;
    QFutureWatcher<search::FindPairResult> find_pair_watcher_;
    QProgressDialog find_pair_progress_;
    GlossaryDialog glossary_;
    CriteriaDialog criteria_dialog_;
    SystemDefaultsDialog defaults_dialog_;
    AboutDialog about_;
};

}  // namespace primerpicker::gui
