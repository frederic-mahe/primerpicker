#pragma once

#include "core/analysis/legacy_verdict.hpp"

#include <QAbstractTableModel>
#include <QFormLayout>
#include <QLabel>
#include <QModelIndex>
#include <QObject>
#include <QPushButton>
#include <QString>
#include <QTableView>
#include <QVariant>
#include <QVBoxLayout>
#include <QWidget>

#include <vector>

namespace primerpicker::gui {

// A row of the criteria table: criterion, setting, actual value and
// whether the primer meets it
struct CriterionRow {
    QString criterion;
    QString setting;
    QString actual;
    bool meets = false;
};

// The criteria table of the Primer Analysis screen, as a Qt model
class CriteriaModel final : public QAbstractTableModel {
    Q_OBJECT

public:
    using QAbstractTableModel::QAbstractTableModel;

    void set_rows(std::vector<CriterionRow> rows);

    [[nodiscard]] auto rowCount(QModelIndex const& parent) const -> int override;
    [[nodiscard]] auto columnCount(QModelIndex const& parent) const -> int override;
    [[nodiscard]] auto data(QModelIndex const& index, int role) const -> QVariant override;
    [[nodiscard]] auto headerData(int section, Qt::Orientation orientation, int role) const -> QVariant override;

private:
    std::vector<CriterionRow> rows_;
};

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
    CriteriaModel criteria_model_;
    QTableView criteria_table_;
    QLabel hairpin_;
    QLabel hairpin_drawing_;
    QLabel verdict_;
    QPushButton print_;
};

}  // namespace primerpicker::gui
