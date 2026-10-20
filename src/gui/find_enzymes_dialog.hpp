#pragma once

#include "core/enzymes/site_search.hpp"

#include <QCheckBox>
#include <QCompleter>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QLabel>
#include <QLineEdit>
#include <QObject>
#include <QString>
#include <QStringListModel>
#include <QVBoxLayout>
#include <QWidget>

namespace primerpicker::gui {

// Molecule-Ops -> Find Enzymes: the enzyme (a REBASE name; old names
// such as Asp718 or Sau3A are accepted) and the kind of sites, true or
// almost (one mismatch), as PRIMER.exe's True Sites and Almost Sites
class FindEnzymesDialog final : public QDialog {
    Q_OBJECT

public:
    explicit FindEnzymesDialog(QWidget& parent);

    // opens the dialog (window-modal, returns at once), with the last
    // enzyme and kind of sites
    void start();

    [[nodiscard]] auto enzyme() const -> QString;
    [[nodiscard]] auto mismatches() const -> enzymes::Mismatches;

    // keeps the dialog open with a message
    void refuse(QString const& message);

Q_SIGNALS:
    // OK was chosen: the main window searches, then accepts the dialog
    // or refuses the enzyme
    void enzyme_entered();

private:
    // members are destroyed in reverse order: the widgets before the
    // layouts that hold them
    QVBoxLayout layout_;
    QFormLayout form_;
    QStringListModel names_;
    QCompleter completer_;
    QLineEdit enzyme_;
    QCheckBox almost_;
    QLabel message_;
    QDialogButtonBox buttons_;
};

}  // namespace primerpicker::gui
