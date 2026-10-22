#pragma once

#include "core/sequence/strand.hpp"

#include <QCheckBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QLabel>
#include <QLineEdit>
#include <QObject>
#include <QRegularExpressionValidator>
#include <QString>
#include <QVBoxLayout>
#include <QWidget>

#include <cstddef>
#include <optional>

namespace primerpicker::gui {

// The search of Primer/Create -> Match Active Primer: primers B from
// target - start_offset to target + end_offset, on the complement strand
// (the form's default) or on the top strand
struct MatchRequest {
    std::size_t target = 0;
    std::size_t start_offset = 0;
    std::size_t end_offset = 0;
    sequence::Strand strand = sequence::Strand::complement;  // of primers B
};

// Primer/Create -> Match Active Primer: the "Define Matching Primer"
// form of PRIMER.exe (the target, and the search range from the
// settings)
class MatchDialog final : public QDialog {
    Q_OBJECT

public:
    explicit MatchDialog(QWidget& parent);

    // opens the dialog with the settings' range and primers B on
    // STRAND_B (window-modal, returns at once)
    void start(std::size_t start_offset, std::size_t end_offset, sequence::Strand strand_b);

    // the typed search, if every field holds a number (target: 1 or more)
    [[nodiscard]] auto request() const -> std::optional<MatchRequest>;

    // keeps the dialog open with a message
    void refuse(QString const& message);

Q_SIGNALS:
    // OK was chosen: the main window searches, then accepts the dialog
    // or refuses the search
    void search_entered();

private:
    // members are destroyed in reverse order: the widgets before the
    // layouts that hold them
    QVBoxLayout layout_;
    QFormLayout form_;
    QRegularExpressionValidator digits_;
    QLineEdit target_;
    QLineEdit start_offset_;
    QLineEdit end_offset_;
    QCheckBox complement_;
    QLabel message_;
    QDialogButtonBox buttons_;
};

}  // namespace primerpicker::gui
