#pragma once

#include "core/sequence/nucleotide.hpp"

#include <QDialog>
#include <QDialogButtonBox>
#include <QLineEdit>
#include <QObject>
#include <QString>
#include <QStringView>
#include <QtTypes>
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

// What a sequence dialog is for: its window title, the prefix of its
// widgets' objectNames (NAME_dialog, NAME_field, NAME_buttons) and the
// fewest bases it accepts
struct SequenceDialogKind {
    QStringView title;
    QStringView name;
    qsizetype min_length;
};

// Primer/Create -> Direct Entry: the sequence of a new primer (at least
// two bases; the analyses need them)
constexpr auto direct_entry_kind =
    SequenceDialogKind{.title = u"Direct Entry", .name = u"direct_entry", .min_length = 2};
// Modify -> 5'-extension: the 5' extension of the active primer (none
// when empty)
constexpr auto extension_kind = SequenceDialogKind{.title = u"5' Extension", .name = u"extension", .min_length = 0};

// A field for a sequence that behaves as PRIMER.exe's Direct Entry
// field (at most 50 bases, A, C, G and T only)
class SequenceDialog final : public QDialog {
    Q_OBJECT

public:
    SequenceDialog(QWidget& parent, SequenceDialogKind const& kind);

    // opens the dialog with TEXT in the field (window-modal, returns at
    // once)
    void start(QString const& text = {});

    // the typed sequence, 5' to 3'
    [[nodiscard]] auto sequence() const -> std::vector<sequence::Nucleotide>;

private:
    void update_ok_button();

    qsizetype min_length_;

    // members are destroyed in reverse order: the widgets before the
    // layout that holds them
    QVBoxLayout layout_;
    NucleotideValidator validator_;
    QLineEdit field_;
    QDialogButtonBox buttons_;
};

}  // namespace primerpicker::gui
