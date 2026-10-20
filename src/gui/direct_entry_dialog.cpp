#include "gui/direct_entry_dialog.hpp"

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

namespace {

// the analyses need at least two bases
constexpr auto min_primer_length = qsizetype{2};

}  // namespace

auto NucleotideValidator::validate(QString& input, [[maybe_unused]] int& position) const -> State {
    input = input.toUpper();
    auto const inserted = std::ranges::all_of(
        input, [](QChar const letter) { return sequence::direct_entry_inserts(letter.toLatin1()); });
    return inserted ? Acceptable : Invalid;
}

DirectEntryDialog::DirectEntryDialog(QWidget& parent)
    : QDialog{&parent},
      validator_{this},
      buttons_{QDialogButtonBox::Ok | QDialogButtonBox::Cancel} {
    setObjectName(u"direct_entry_dialog"_s);
    setWindowTitle(u"Direct Entry"_s);
    field_.setObjectName(u"direct_entry_field"_s);
    field_.setMaxLength(static_cast<int>(sequence::direct_entry_max_length));
    field_.setValidator(&validator_);
    buttons_.setObjectName(u"direct_entry_buttons"_s);
    layout_.addWidget(&field_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    QObject::connect(&field_, &QLineEdit::textChanged, this, &DirectEntryDialog::update_ok_button);
    QObject::connect(&buttons_, &QDialogButtonBox::accepted, this, &QDialog::accept);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
    update_ok_button();
}

void DirectEntryDialog::start() {
    field_.clear();
    field_.setFocus();
    open();
}

auto DirectEntryDialog::primer() const -> std::vector<sequence::Nucleotide> {
    auto primer = std::vector<sequence::Nucleotide>{};
    for (auto const letter : field_.text()) {
        auto const nucleotide = sequence::to_nucleotide(letter.toLatin1());
        if (nucleotide.has_value()) {  // always: the validator inserts nucleotides only
            primer.push_back(*nucleotide);
        }
    }
    return primer;
}

void DirectEntryDialog::update_ok_button() {
    buttons_.button(QDialogButtonBox::Ok)->setEnabled(field_.text().size() >= min_primer_length);
}

}  // namespace primerpicker::gui
