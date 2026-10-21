#include "gui/multi_primer_dialog.hpp"

#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QList>
#include <QObject>
#include <QString>
#include <QWidget>

#include <cstddef>
#include <optional>
#include <utility>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

constexpr auto field_names = QStringView{u"abcd"};
constexpr auto empty_choice = 0;  // the first item of a field: empty

}  // namespace

MultiPrimerDialog::MultiPrimerDialog(QWidget& parent)
    : QDialog{&parent}, buttons_{QDialogButtonBox::Ok | QDialogButtonBox::Cancel} {
    setObjectName(u"multi_primer_dialog"_s);
    setWindowTitle(u"Multi Primer"_s);
    buttons_.setObjectName(u"multi_primer_buttons"_s);
    for (auto index = std::size_t{0}; index < fields_.size(); ++index) {
        auto const name = field_names.at(static_cast<qsizetype>(index));
        fields_.at(index).setObjectName(u"multi_%1_box"_s.arg(name));
        form_.addRow(u"Primer &%1:"_s.arg(name.toUpper()), &fields_.at(index));
    }
    layout_.addLayout(&form_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    QObject::connect(&buttons_, &QDialogButtonBox::accepted, this, &QDialog::accept);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
}

void MultiPrimerDialog::start(QList<QString> const& primers) {
    for (auto index = std::size_t{0}; index < fields_.size(); ++index) {
        auto& field = fields_.at(index);
        field.clear();
        field.addItem(u"(empty)"_s);
        field.addItems(primers);
        // the first primers of the list, one per field
        field.setCurrentIndex(std::cmp_less(index, primers.size()) ? static_cast<int>(index) + 1 : empty_choice);
    }
    fields_.front().setFocus();
    open();
}

auto MultiPrimerDialog::selection() const -> MultiSelection {
    auto selection = MultiSelection{};
    for (auto index = std::size_t{0}; index < fields_.size(); ++index) {
        auto const choice = fields_.at(index).currentIndex();
        if (choice > empty_choice) {
            selection.at(index) = static_cast<std::size_t>(choice - 1);
        }
    }
    return selection;
}

}  // namespace primerpicker::gui
