#include "gui/criteria_dialog.hpp"

#include "core/analysis/legacy_verdict.hpp"

#include <QDialog>
#include <QDialogButtonBox>
#include <QObject>
#include <QString>
#include <QWidget>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// PRIMER.exe's message for a refused value (from the strings of
// examples/PRIMER.exe)
constexpr auto illegal_value_message = u"Illegal value for defaults";

}  // namespace

CriteriaDialog::CriteriaDialog(QWidget& parent)
    : QDialog{&parent},
      fields_{u"criteria_dialog"_s},
      buttons_{QDialogButtonBox::Ok | QDialogButtonBox::Cancel} {
    setObjectName(u"criteria_dialog"_s);
    setWindowTitle(u"Criteria of the Active Primer"_s);
    message_.setObjectName(u"criteria_dialog_message"_s);
    buttons_.setObjectName(u"criteria_dialog_buttons"_s);
    layout_.addWidget(&fields_);
    layout_.addWidget(&message_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    QObject::connect(&buttons_, &QDialogButtonBox::accepted, this, &CriteriaDialog::check_and_accept);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
}

void CriteriaDialog::start(analysis::LegacyCriteria const& criteria) {
    fields_.set_criteria(criteria);
    message_.clear();
    open();
}

auto CriteriaDialog::criteria() const -> analysis::LegacyCriteria {
    return fields_.criteria();
}

void CriteriaDialog::check_and_accept() {
    if (not fields_.valid()) {
        message_.setText(QString::fromUtf16(illegal_value_message));
        return;
    }
    accept();
    Q_EMIT criteria_entered();
}

}  // namespace primerpicker::gui
