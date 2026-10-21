#include "gui/match_dialog.hpp"

#include "core/sequence/strand.hpp"

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

MatchDialog::MatchDialog(QWidget& parent)
    : QDialog{&parent},
      digits_{QRegularExpression{u"[0-9]{1,18}"_s}},
      complement_{u"Primer B on the &complement strand"_s},
      buttons_{QDialogButtonBox::Ok | QDialogButtonBox::Cancel} {
    setObjectName(u"match_dialog"_s);
    setWindowTitle(u"Match Active Primer"_s);
    target_.setObjectName(u"match_target_field"_s);
    start_offset_.setObjectName(u"match_start_field"_s);
    end_offset_.setObjectName(u"match_end_field"_s);
    message_.setObjectName(u"match_message_label"_s);
    complement_.setObjectName(u"match_complement_box"_s);
    complement_.setChecked(true);
    for (auto* const field : std::array{&target_, &start_offset_, &end_offset_}) {
        field->setValidator(&digits_);
        QObject::connect(field, &QLineEdit::textChanged, &message_, &QLabel::clear);
    }
    form_.addRow(u"Primer B, &target bp position:"_s, &target_);
    form_.addRow(u"    &Start: Target -"_s, &start_offset_);
    form_.addRow(u"    &End: Target +"_s, &end_offset_);
    form_.addRow(&complement_);
    layout_.addLayout(&form_);
    layout_.addWidget(&message_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    QObject::connect(&buttons_, &QDialogButtonBox::accepted, this, &MatchDialog::search_entered);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
}

void MatchDialog::start(std::size_t const start_offset, std::size_t const end_offset) {
    start_offset_.setText(QString::number(start_offset));
    end_offset_.setText(QString::number(end_offset));
    message_.clear();
    target_.selectAll();
    target_.setFocus();
    open();
}

auto MatchDialog::request() const -> std::optional<MatchRequest> {
    auto const target = number(target_);
    auto const start_offset = number(start_offset_);
    auto const end_offset = number(end_offset_);
    if (not target or not start_offset or not end_offset or *target == 0) {
        return std::nullopt;
    }
    return MatchRequest{.target = *target,
                        .start_offset = *start_offset,
                        .end_offset = *end_offset,
                        .strand = complement_.isChecked() ? sequence::Strand::complement : sequence::Strand::top};
}

void MatchDialog::refuse(QString const& message) {
    message_.setText(message);
    target_.setFocus();
}

}  // namespace primerpicker::gui
