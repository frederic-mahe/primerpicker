#pragma once

#include "core/sequence/strand.hpp"

#include <QCheckBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QLabel>
#include <QLineEdit>
#include <QObject>
#include <QRegularExpressionValidator>
#include <QString>
#include <QVBoxLayout>
#include <QWidget>

#include <cstddef>
#include <optional>

namespace primerpicker::gui {

// Primer/Create -> From Molecule: the base position of a new primer on
// the loaded molecule, and its strand, as PRIMER.exe asks for them
// ("Base position for primer:", "Complement strand (Y/N):"). The
// length is the System Defaults' primer length.
class FromMoleculeDialog final : public QDialog {
    Q_OBJECT

public:
    explicit FromMoleculeDialog(QWidget& parent);

    // opens the dialog for MOLECULE_NAME (window-modal, returns at once)
    void start(QString const& molecule_name);

    // the typed position (1 or more), if any
    [[nodiscard]] auto position() const -> std::optional<std::size_t>;
    [[nodiscard]] auto strand() const -> sequence::Strand;

    // keeps the dialog open with PRIMER.exe's message
    void refuse_position();

Q_SIGNALS:
    // OK was chosen: the main window makes the primer, then accepts the
    // dialog or refuses the position
    void position_entered();

private:
    // members are destroyed in reverse order: the widgets before the
    // layouts that hold them
    QVBoxLayout layout_;
    QFormLayout form_;
    QRegularExpressionValidator digits_;
    QLabel molecule_;
    QLineEdit position_;
    QCheckBox complement_;
    QLabel message_;
    QDialogButtonBox buttons_;
};

}  // namespace primerpicker::gui
