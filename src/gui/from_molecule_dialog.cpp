#include "gui/from_molecule_dialog.hpp"

#include "core/sequence/strand.hpp"

#include <QDialog>
#include <QDialogButtonBox>
#include <QObject>
#include <QPushButton>
#include <QRegularExpression>
#include <QString>
#include <QWidget>

#include <cstddef>
#include <optional>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

FromMoleculeDialog::FromMoleculeDialog(QWidget& parent)
    : QDialog{&parent},
      digits_{QRegularExpression{u"[0-9]{1,18}"_s}},
      complement_{u"&Complement strand"_s},
      buttons_{QDialogButtonBox::Ok | QDialogButtonBox::Cancel} {
    setObjectName(u"from_molecule_dialog"_s);
    setWindowTitle(u"From Molecule"_s);
    position_.setObjectName(u"position_field"_s);
    position_.setValidator(&digits_);
    complement_.setObjectName(u"complement_box"_s);
    message_.setObjectName(u"position_message"_s);
    form_.addRow(u"From molecule:"_s, &molecule_);
    form_.addRow(u"&Base position for primer:"_s, &position_);
    form_.addRow(&complement_);
    layout_.addLayout(&form_);
    layout_.addWidget(&message_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    QObject::connect(&buttons_, &QDialogButtonBox::accepted, this, &FromMoleculeDialog::position_entered);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
    QObject::connect(&position_, &QLineEdit::textChanged, &message_, &QLabel::clear);
}

void FromMoleculeDialog::start(QString const& molecule_name) {
    molecule_.setText(molecule_name);
    position_.clear();
    message_.clear();
    position_.setFocus();
    open();
}

auto FromMoleculeDialog::position() const -> std::optional<std::size_t> {
    auto converted = false;
    auto const value = position_.text().toULongLong(&converted);
    if (not converted or value == 0) {
        return std::nullopt;
    }
    return std::size_t{value};
}

auto FromMoleculeDialog::strand() const -> sequence::Strand {
    return complement_.isChecked() ? sequence::Strand::complement : sequence::Strand::top;
}

void FromMoleculeDialog::refuse_position() {
    message_.setText(u"Invalid Position"_s);  // PRIMER.exe's message
    position_.setFocus();
}

}  // namespace primerpicker::gui
