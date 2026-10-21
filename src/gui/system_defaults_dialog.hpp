#pragma once

#include "core/io/pd2_cfg.hpp"
#include "gui/criteria_fields.hpp"

#include <QCheckBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QLabel>
#include <QObject>
#include <QPushButton>
#include <QSpinBox>
#include <QString>
#include <QVBoxLayout>
#include <QWidget>

namespace primerpicker::gui {

// Modify -> System Defaults: PRIMER.exe's System Defaults, Criteria and
// Find-Pair weighting factors screens in one dialog. A minimum above
// its maximum is refused with PRIMER.exe's "Illegal value for
// defaults". The dialog also loads and saves settings files (pd2.cfg).
class SystemDefaultsDialog final : public QDialog {
    Q_OBJECT

public:
    explicit SystemDefaultsDialog(QWidget& parent);

    // opens the dialog with DEFAULTS (window-modal, returns at once)
    void start(io::LegacySystemDefaults const& defaults);

    // the values of the fields (bytes of unknown meaning as started)
    [[nodiscard]] auto defaults() const -> io::LegacySystemDefaults;

Q_SIGNALS:
    // OK was chosen with valid values (the dialog is then closed)
    void defaults_entered();
    void load_requested();
    void save_requested();

private:
    void check_and_accept();

    io::LegacySystemDefaults started_;
    // members are destroyed in reverse order: the fields, then the layout
    // that holds them, then the box that holds the layout (a box deletes
    // its layout: the layout must go first)
    QVBoxLayout layout_;
    QGroupBox system_box_;
    QFormLayout system_layout_;
    QSpinBox primer_length_;
    QCheckBox exclude_;
    QCheckBox sort_;
    QSpinBox gc_tolerance_;
    QSpinBox tm_tolerance_;
    QSpinBox match_start_;
    QSpinBox match_end_;
    QSpinBox find_pair_start_;
    QSpinBox find_pair_end_;
    QGroupBox criteria_box_;
    QVBoxLayout criteria_layout_;
    CriteriaFields criteria_;
    QGroupBox weights_box_;
    QFormLayout weights_layout_;
    QSpinBox weight_gc_;
    QSpinBox weight_tm_;
    QSpinBox weight_3prime_;
    QSpinBox weight_any_;
    QSpinBox weight_false_priming_;
    QLabel message_;
    QDialogButtonBox buttons_;
    QPushButton load_;
    QPushButton save_;
};

}  // namespace primerpicker::gui
