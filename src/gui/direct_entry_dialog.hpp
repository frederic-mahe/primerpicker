#pragma once

#include "core/sequence/nucleotide.hpp"

#include <QDialog>
#include <QDialogButtonBox>
#include <QLineEdit>
#include <QObject>
#include <QString>
#include <QValidator>
#include <QVBoxLayout>
#include <QWidget>

#include <vector>

namespace primerpicker::gui {

// The letters PRIMER.exe's Direct Entry field inserts: A, C, G and T
// (typed lower-case letters are shown upper-case); anything else is
// not inserted
class NucleotideValidator final : public QValidator {
    Q_OBJECT

public:
    using QValidator::QValidator;

    auto validate(QString& input, int& position) const -> State override;
};

// Primer/Create -> Direct Entry: types in the sequence of a new primer
// (at most 50 bases, at least two)
class DirectEntryDialog final : public QDialog {
    Q_OBJECT

public:
    explicit DirectEntryDialog(QWidget& parent);

    // opens the dialog with an empty field (window-modal, returns at once)
    void start();

    // the typed primer, 5' to 3'
    [[nodiscard]] auto primer() const -> std::vector<sequence::Nucleotide>;

private:
    void update_ok_button();

    // members are destroyed in reverse order: the widgets before the
    // layout that holds them
    QVBoxLayout layout_;
    NucleotideValidator validator_;
    QLineEdit field_;
    QDialogButtonBox buttons_;
};

}  // namespace primerpicker::gui
