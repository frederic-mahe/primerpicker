#pragma once

#include "core/analysis/legacy_verdict.hpp"
#include "gui/text_table_model.hpp"

#include <QLabel>
#include <QObject>
#include <QTableView>
#include <QVBoxLayout>
#include <QWidget>

namespace primerpicker::gui {

// The Primer Pair Analysis screen of PRIMER.exe: primer balance (length,
// % GC, Tm and their differences), cautions (runs, hairpins) and the
// A:A, B:B and A:B dimers, with PRIMER.exe's "!" flags
class PairView final : public QWidget {
    Q_OBJECT

public:
    PairView();

    void show_pair(analysis::PrimerPair const& primers, analysis::LegacyCriteria const& criteria);

private:
    // members are destroyed in reverse order: the widgets before the
    // layout that holds them
    QVBoxLayout layout_;
    QLabel primer_a_;
    QLabel primer_b_;
    TextTableModel model_;
    QTableView table_;
};

}  // namespace primerpicker::gui
