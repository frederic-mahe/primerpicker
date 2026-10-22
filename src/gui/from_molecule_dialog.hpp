#pragma once

#include "core/sequence/nucleotide.hpp"
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
#include <vector>

namespace primerpicker::gui {

// Primer/Create -> From Molecule: the base position of a new primer on
// the loaded molecule, its strand, its length (5' extension included)
// and its 5' extension, as PRIMER.exe asks for them ("Base position for
// primer:", "Complement strand (Y/N):", "Primer length (incl 5' ext):",
// "5' extension to be added:"). The length starts at the System
// Defaults' primer length.
class FromMoleculeDialog final : public QDialog {
    Q_OBJECT

public:
    explicit FromMoleculeDialog(QWidget& parent);

    // opens the dialog for MOLECULE_NAME, with LENGTH in the length field
    // (window-modal, returns at once)
    void start(QString const& molecule_name, std::size_t length);

    // the typed position (1 or more), if any
    [[nodiscard]] auto position() const -> std::optional<std::size_t>;
    [[nodiscard]] auto strand() const -> sequence::Strand;
    // the typed length, extension included (1 to 99), if any
    [[nodiscard]] auto length() const -> std::optional<std::size_t>;
    [[nodiscard]] auto extension() const -> std::vector<sequence::Nucleotide>;

    // keeps the dialog open with PRIMER.exe's message
    void refuse_position();
    // keeps the dialog open with MESSAGE
    void refuse(QString const& message);

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
    QRegularExpressionValidator length_digits_;
    QRegularExpressionValidator bases_;
    QLabel molecule_;
    QLineEdit position_;
    QCheckBox complement_;
    QLineEdit length_;
    QLineEdit extension_;
    QLabel message_;
    QDialogButtonBox buttons_;
};

}  // namespace primerpicker::gui
