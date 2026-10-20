#include "gui/sequence_dialog.hpp"

#include "core/sequence/legacy_direct_entry.hpp"
#include "core/sequence/nucleotide.hpp"

#include <QDialog>
#include <QDialogButtonBox>
#include <QLabel>
#include <QLineEdit>
#include <QObject>
#include <QPushButton>
#include <QString>
#include <QValidator>
#include <QWidget>

#include <algorithm>
#include <vector>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

auto NucleotideValidator::validate(QString& input, [[maybe_unused]] int& position) const -> State {
    input = input.toUpper();
    auto const inserted = std::ranges::all_of(
        input, [](QChar const letter) { return sequence::direct_entry_inserts(letter.toLatin1()); });
    return inserted ? Acceptable : Invalid;
}

SequenceDialog::SequenceDialog(QWidget& parent, SequenceDialogKind const& kind)
    : QDialog{&parent},
      min_length_{kind.min_length},
      validator_{this},
      buttons_{QDialogButtonBox::Ok | QDialogButtonBox::Cancel} {
    auto const name = kind.name.toString();
    setObjectName(name + u"_dialog"_s);
    setWindowTitle(kind.title.toString());
    field_.setObjectName(name + u"_field"_s);
    field_.setMaxLength(static_cast<int>(sequence::direct_entry_max_length));
    field_.setValidator(&validator_);
    buttons_.setObjectName(name + u"_buttons"_s);
    layout_.addWidget(&field_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    QObject::connect(&field_, &QLineEdit::textChanged, this, &SequenceDialog::update_ok_button);
    QObject::connect(&buttons_, &QDialogButtonBox::accepted, this, &QDialog::accept);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
    update_ok_button();
}

void SequenceDialog::start(QString const& text) {
    field_.setText(text);
    field_.setFocus();
    open();
}

auto SequenceDialog::sequence() const -> std::vector<sequence::Nucleotide> {
    auto nucleotides = std::vector<sequence::Nucleotide>{};
    for (auto const letter : field_.text()) {
        auto const nucleotide = sequence::to_nucleotide(letter.toLatin1());
        if (nucleotide.has_value()) {  // always: the validator inserts nucleotides only
            nucleotides.push_back(*nucleotide);
        }
    }
    return nucleotides;
}

void SequenceDialog::update_ok_button() {
    buttons_.button(QDialogButtonBox::Ok)->setEnabled(field_.text().size() >= min_length_);
}

}  // namespace primerpicker::gui
