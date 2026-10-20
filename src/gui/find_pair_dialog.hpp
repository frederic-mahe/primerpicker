#pragma once

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

// The search of Primer/Create -> Find Pair: primer A (top strand) from
// target_a - start_offset_a to target_a, primer B (complement strand)
// from target_b to target_b + end_offset_b
struct FindPairRequest {
    std::size_t target_a = 0;
    std::size_t start_offset_a = 0;
    std::size_t target_b = 0;
    std::size_t end_offset_b = 0;
};

// Primer/Create -> Find Pair: the "Define Primer Pair Constraints" form
// of PRIMER.exe (targets, and the search range, from the settings)
class FindPairDialog final : public QDialog {
    Q_OBJECT

public:
    explicit FindPairDialog(QWidget& parent);

    // opens the dialog with the settings' offsets (window-modal, returns
    // at once)
    void start(std::size_t start_offset_a, std::size_t end_offset_b);

    // the typed search, if every field holds a number (targets: 1 or
    // more)
    [[nodiscard]] auto request() const -> std::optional<FindPairRequest>;

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
    QLineEdit target_a_;
    QLineEdit start_offset_;
    QLineEdit target_b_;
    QLineEdit end_offset_;
    QLabel message_;
    QDialogButtonBox buttons_;
};

}  // namespace primerpicker::gui
