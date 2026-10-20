#pragma once

#include "core/io/pd2_cfg.hpp"
#include "core/sequence/nucleotide.hpp"
#include "gui/analysis_view.hpp"
#include "gui/direct_entry_dialog.hpp"

#include <QLabel>
#include <QMainWindow>
#include <QObject>
#include <QStackedWidget>

#include <cstddef>
#include <optional>
#include <vector>

namespace primerpicker::gui {

// The primer the analyses work on: the whole sequence 5' to 3', whose
// first extension_length nucleotides are the 5' extension
struct ActivePrimer {
    std::vector<sequence::Nucleotide> sequence;
    std::size_t extension_length = 0;
};

// The main window: Primer Designer's main menu as a menu bar. Items
// that are not implemented yet are disabled.
class MainWindow final : public QMainWindow {
    Q_OBJECT

public:
    MainWindow();

private:
    void add_legacy_menus();
    void use_direct_entry();
    void analyze_active_primer();

    io::LegacySystemDefaults defaults_;
    std::optional<ActivePrimer> active_primer_;
    // members are destroyed in reverse order: the pages before the
    // stack that holds them
    QStackedWidget pages_;
    QLabel no_primer_;
    AnalysisView analysis_;
    DirectEntryDialog direct_entry_;
};

}  // namespace primerpicker::gui
