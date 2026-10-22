#include "gui/from_molecule_dialog.hpp"

#include "core/sequence/nucleotide.hpp"
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
#include <vector>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

FromMoleculeDialog::FromMoleculeDialog(QWidget& parent)
    : QDialog{&parent},
      digits_{QRegularExpression{u"[0-9]{1,18}"_s}},
      // PRIMER.exe's field holds 2 digits (tests/legacy/README.md, "The
      // picking screen")
      length_digits_{QRegularExpression{u"[0-9]{1,2}"_s}},
      bases_{QRegularExpression{u"[ACGTacgt]*"_s}},
      complement_{u"&Complement strand"_s},
      buttons_{QDialogButtonBox::Ok | QDialogButtonBox::Cancel} {
    setObjectName(u"from_molecule_dialog"_s);
    setWindowTitle(u"From Molecule"_s);
    position_.setObjectName(u"position_field"_s);
    position_.setValidator(&digits_);
    complement_.setObjectName(u"complement_box"_s);
    length_.setObjectName(u"length_field"_s);
    length_.setValidator(&length_digits_);
    extension_.setObjectName(u"extension_field_form"_s);
    extension_.setValidator(&bases_);
    message_.setObjectName(u"position_message"_s);
    form_.addRow(u"From molecule:"_s, &molecule_);
    form_.addRow(u"&Base position for primer:"_s, &position_);
    form_.addRow(&complement_);
    form_.addRow(u"Primer &length (incl 5' ext):"_s, &length_);
    form_.addRow(u"5' e&xtension to be added:"_s, &extension_);
    layout_.addLayout(&form_);
    layout_.addWidget(&message_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    QObject::connect(&buttons_, &QDialogButtonBox::accepted, this, &FromMoleculeDialog::position_entered);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
    for (auto* const field : {&position_, &length_, &extension_}) {
        QObject::connect(field, &QLineEdit::textChanged, &message_, &QLabel::clear);
    }
}

void FromMoleculeDialog::start(QString const& molecule_name, std::size_t const length) {
    molecule_.setText(molecule_name);
    position_.clear();
    length_.setText(QString::number(length));
    extension_.clear();
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

auto FromMoleculeDialog::length() const -> std::optional<std::size_t> {
    auto converted = false;
    auto const value = length_.text().toULongLong(&converted);
    if (not converted or value == 0) {
        return std::nullopt;
    }
    return std::size_t{value};
}

auto FromMoleculeDialog::extension() const -> std::vector<sequence::Nucleotide> {
    // the validator lets letters A, C, G and T only
    return sequence::to_nucleotides(extension_.text().toUpper().toStdString()).value_or(std::vector<sequence::Nucleotide>{});
}

void FromMoleculeDialog::refuse(QString const& message) {
    message_.setText(message);
    length_.setFocus();
}

void FromMoleculeDialog::refuse_position() {
    message_.setText(u"Invalid Position"_s);  // PRIMER.exe's message
    position_.setFocus();
}

}  // namespace primerpicker::gui
