#pragma once

#include "core/analysis/legacy_verdict.hpp"
#include "gui/criteria_fields.hpp"

#include <QDialog>
#include <QDialogButtonBox>
#include <QLabel>
#include <QObject>
#include <QString>
#include <QVBoxLayout>
#include <QWidget>

namespace primerpicker::gui {

// Modify -> Active Primer: Criteria: the criteria of the active primer
// (PRIMER.exe's Criteria form); a minimum above its maximum is refused
// with PRIMER.exe's "Illegal value for defaults"
class CriteriaDialog final : public QDialog {
    Q_OBJECT

public:
    explicit CriteriaDialog(QWidget& parent);

    // opens the dialog with CRITERIA (window-modal, returns at once)
    void start(analysis::LegacyCriteria const& criteria);

    [[nodiscard]] auto criteria() const -> analysis::LegacyCriteria;

Q_SIGNALS:
    // OK was chosen with valid values (the dialog is then closed)
    void criteria_entered();

private:
    void check_and_accept();

    // members are destroyed in reverse order: the widgets before the
    // layout that holds them
    QVBoxLayout layout_;
    CriteriaFields fields_;
    QLabel message_;
    QDialogButtonBox buttons_;
};

}  // namespace primerpicker::gui
