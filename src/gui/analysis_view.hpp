#pragma once

#include "core/analysis/legacy_verdict.hpp"
#include "gui/text_table_model.hpp"

#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QObject>
#include <QPushButton>
#include <QTableView>
#include <QVBoxLayout>
#include <QWidget>

namespace primerpicker::gui {

// Analyze -> Active Primer: the Primer Analysis screen of PRIMER.exe
// (sequence, length, 5' extension, composition, the criteria table, the
// most stable hairpin and the overall verdict)
class AnalysisView final : public QWidget {
    Q_OBJECT

public:
    AnalysisView();

    void show_primer(analysis::LegacyPrimer const& primer, analysis::LegacyCriteria const& criteria);

Q_SIGNALS:
    // Print to Disk was chosen (P, as on PRIMER.exe's screen)
    void print_requested();
    // Uniqueness (U) or False Priming (F) was chosen
    void uniqueness_requested();
    void false_priming_requested();

private:
    // members are destroyed in reverse order: the widgets before the
    // layouts that hold them
    QVBoxLayout layout_;
    QFormLayout primer_layout_;
    QFormLayout result_layout_;
    QLabel sequence_;
    QLabel length_;
    QLabel extension_;
    QLabel composition_;
    TextTableModel criteria_model_;
    QTableView criteria_table_;
    QLabel hairpin_;
    QLabel hairpin_drawing_;
    QLabel verdict_;
    QHBoxLayout buttons_layout_;
    QPushButton print_;
    QPushButton uniqueness_;
    QPushButton false_priming_;
};

}  // namespace primerpicker::gui
