#pragma once

#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QList>
#include <QObject>
#include <QString>
#include <QVBoxLayout>
#include <QWidget>

#include <cstddef>

namespace primerpicker::gui {

// Two primers of the Primer List, by their place in it (0 for the first)
struct PairSelection {
    std::size_t first = 0;
    std::size_t second = 0;
};

// Analyze -> Primer Pair: primers A and B, chosen in the Primer List (as
// PRIMER.exe asks for their numbers)
class PairDialog final : public QDialog {
    Q_OBJECT

public:
    explicit PairDialog(QWidget& parent);

    // opens the dialog with the primers of the list, one text each
    // (window-modal, returns at once)
    void start(QList<QString> const& primers);

    [[nodiscard]] auto selection() const -> PairSelection;

private:
    // members are destroyed in reverse order: the widgets before the
    // layouts that hold them
    QVBoxLayout layout_;
    QFormLayout form_;
    QComboBox primer_a_;
    QComboBox primer_b_;
    QDialogButtonBox buttons_;
};

}  // namespace primerpicker::gui
