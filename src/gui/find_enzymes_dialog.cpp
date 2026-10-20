#include "gui/find_enzymes_dialog.hpp"

#include "core/enzymes/rebase_enzymes.hpp"
#include "core/enzymes/site_search.hpp"

#include <QCompleter>
#include <QDialog>
#include <QDialogButtonBox>
#include <QLabel>
#include <QList>
#include <QObject>
#include <QString>
#include <QWidget>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// the commercially available enzymes, for completion
auto commercial_names() -> QList<QString> {
    auto names = QList<QString>{};
    for (auto const& enzyme : enzymes::rebase_enzymes) {
        names.append(QString::fromLatin1(enzyme.name));
    }
    return names;
}

}  // namespace

FindEnzymesDialog::FindEnzymesDialog(QWidget& parent)
    : QDialog{&parent},
      names_{commercial_names()},
      almost_{u"&Almost sites (one mismatch)"_s},
      buttons_{QDialogButtonBox::Ok | QDialogButtonBox::Cancel} {
    setObjectName(u"find_enzymes_dialog"_s);
    setWindowTitle(u"Find Enzymes"_s);
    enzyme_.setObjectName(u"enzyme_field"_s);
    almost_.setObjectName(u"almost_box"_s);
    message_.setObjectName(u"enzyme_message"_s);
    completer_.setModel(&names_);
    completer_.setCaseSensitivity(Qt::CaseInsensitive);
    enzyme_.setCompleter(&completer_);
    form_.addRow(u"&Enzyme:"_s, &enzyme_);
    form_.addRow(&almost_);
    layout_.addLayout(&form_);
    layout_.addWidget(&message_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    QObject::connect(&buttons_, &QDialogButtonBox::accepted, this, &FindEnzymesDialog::enzyme_entered);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
    QObject::connect(&enzyme_, &QLineEdit::textChanged, &message_, &QLabel::clear);
}

void FindEnzymesDialog::start() {
    message_.clear();
    enzyme_.selectAll();
    enzyme_.setFocus();
    open();
}

auto FindEnzymesDialog::enzyme() const -> QString {
    return enzyme_.text().trimmed();
}

auto FindEnzymesDialog::mismatches() const -> enzymes::Mismatches {
    return almost_.isChecked() ? enzymes::Mismatches::at_most_one : enzymes::Mismatches::none;
}

void FindEnzymesDialog::refuse(QString const& message) {
    message_.setText(message);
    enzyme_.setFocus();
}

}  // namespace primerpicker::gui
