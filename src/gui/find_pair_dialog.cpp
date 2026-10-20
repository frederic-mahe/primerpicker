#include "gui/find_pair_dialog.hpp"

#include <QDialog>
#include <QDialogButtonBox>
#include <QLabel>
#include <QLineEdit>
#include <QObject>
#include <QRegularExpression>
#include <QString>
#include <QWidget>

#include <array>
#include <cstddef>
#include <optional>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// the number in FIELD, if any
auto number(QLineEdit const& field) -> std::optional<std::size_t> {
    auto converted = false;
    auto const value = field.text().toULongLong(&converted);
    return converted ? std::optional{std::size_t{value}} : std::nullopt;
}

}  // namespace

FindPairDialog::FindPairDialog(QWidget& parent)
    : QDialog{&parent},
      digits_{QRegularExpression{u"[0-9]{1,18}"_s}},
      buttons_{QDialogButtonBox::Ok | QDialogButtonBox::Cancel} {
    setObjectName(u"find_pair_dialog"_s);
    setWindowTitle(u"Find Pair"_s);
    target_a_.setObjectName(u"target_a_field"_s);
    start_offset_.setObjectName(u"start_offset_field"_s);
    target_b_.setObjectName(u"target_b_field"_s);
    end_offset_.setObjectName(u"end_offset_field"_s);
    message_.setObjectName(u"find_pair_message_label"_s);
    for (auto* const field : std::array{&target_a_, &start_offset_, &target_b_, &end_offset_}) {
        field->setValidator(&digits_);
        QObject::connect(field, &QLineEdit::textChanged, &message_, &QLabel::clear);
    }
    form_.addRow(u"Primer A (normal strand), target bp &position:"_s, &target_a_);
    form_.addRow(u"    &Start: Target -"_s, &start_offset_);
    form_.addRow(u"Primer B (compl. strand), target bp p&osition:"_s, &target_b_);
    form_.addRow(u"    &End: Target +"_s, &end_offset_);
    layout_.addLayout(&form_);
    layout_.addWidget(&message_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    QObject::connect(&buttons_, &QDialogButtonBox::accepted, this, &FindPairDialog::search_entered);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
}

void FindPairDialog::start(std::size_t const start_offset_a, std::size_t const end_offset_b) {
    start_offset_.setText(QString::number(start_offset_a));
    end_offset_.setText(QString::number(end_offset_b));
    message_.clear();
    target_a_.selectAll();
    target_a_.setFocus();
    open();
}

auto FindPairDialog::request() const -> std::optional<FindPairRequest> {
    auto const target_a = number(target_a_);
    auto const start_offset = number(start_offset_);
    auto const target_b = number(target_b_);
    auto const end_offset = number(end_offset_);
    if (not target_a or not start_offset or not target_b or not end_offset or *target_a == 0 or *target_b == 0) {
        return std::nullopt;
    }
    return FindPairRequest{
        .target_a = *target_a, .start_offset_a = *start_offset, .target_b = *target_b, .end_offset_b = *end_offset};
}

void FindPairDialog::refuse(QString const& message) {
    message_.setText(message);
    target_a_.setFocus();
}

}  // namespace primerpicker::gui
