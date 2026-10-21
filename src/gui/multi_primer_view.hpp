#pragma once

#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "gui/text_table_model.hpp"

#include <QLabel>
#include <QObject>
#include <QTableView>
#include <QVBoxLayout>
#include <QWidget>

namespace primerpicker::gui {

// The Multi Primer Analysis screen of PRIMER.exe: "Primer Balance and
// Cautions" (a row per filled field) and the two "Primer Dimers"
// matrices, 3' end and Any (self-dimers on the diagonal, pairs above it,
// dots below it)
class MultiPrimerView final : public QWidget {
    Q_OBJECT

public:
    MultiPrimerView();

    void show_primers(analysis::MultiPrimerFields const& fields, analysis::LegacyCriteria const& criteria);

private:
    // members are destroyed in reverse order: the widgets before the
    // layout that holds them
    QVBoxLayout layout_;
    QLabel balance_title_;
    TextTableModel balance_model_;
    QTableView balance_table_;
    QLabel three_prime_title_;
    TextTableModel three_prime_model_;
    QTableView three_prime_table_;
    QLabel any_title_;
    TextTableModel any_model_;
    QTableView any_table_;
};

}  // namespace primerpicker::gui
