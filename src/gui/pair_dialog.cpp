#include "gui/pair_dialog.hpp"

#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QList>
#include <QObject>
#include <QString>
#include <QWidget>

#include <array>
#include <cstddef>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

PairDialog::PairDialog(QWidget& parent)
    : QDialog{&parent}, buttons_{QDialogButtonBox::Ok | QDialogButtonBox::Cancel} {
    setObjectName(u"primer_pair_dialog"_s);
    setWindowTitle(u"Primer Pair"_s);
    primer_a_.setObjectName(u"pair_a_box"_s);
    primer_b_.setObjectName(u"pair_b_box"_s);
    buttons_.setObjectName(u"primer_pair_buttons"_s);
    form_.addRow(u"Primer &A:"_s, &primer_a_);
    form_.addRow(u"Primer &B:"_s, &primer_b_);
    layout_.addLayout(&form_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    QObject::connect(&buttons_, &QDialogButtonBox::accepted, this, &QDialog::accept);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
}

void PairDialog::start(QList<QString> const& primers) {
    for (auto* const box : std::array{&primer_a_, &primer_b_}) {
        box->clear();
        box->addItems(primers);
    }
    // the last two primers saved, as the most likely pair
    constexpr auto pair_size = qsizetype{2};
    if (primers.size() >= pair_size) {
        primer_a_.setCurrentIndex(static_cast<int>(primers.size() - pair_size));
        primer_b_.setCurrentIndex(static_cast<int>(primers.size() - 1));
    }
    primer_a_.setFocus();
    open();
}

auto PairDialog::selection() const -> PairSelection {
    return {.first = static_cast<std::size_t>(primer_a_.currentIndex()),
            .second = static_cast<std::size_t>(primer_b_.currentIndex())};
}

}  // namespace primerpicker::gui
