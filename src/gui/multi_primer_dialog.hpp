#pragma once

#include "core/analysis/legacy_multi_primer.hpp"

#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QList>
#include <QObject>
#include <QString>
#include <QVBoxLayout>
#include <QWidget>

#include <array>
#include <cstddef>
#include <optional>

namespace primerpicker::gui {

// Fields Primer A to Primer D: a primer of the Primer List, by its place
// in it (0 for the first), or empty
using MultiSelection = std::array<std::optional<std::size_t>, analysis::legacy_multi_primer_fields>;

// Analyze -> Multi Primer: up to four primers of the Primer List, in
// fields A to D (as PRIMER.exe asks for their numbers; a field may stay
// empty)
class MultiPrimerDialog final : public QDialog {
    Q_OBJECT

public:
    explicit MultiPrimerDialog(QWidget& parent);

    // opens the dialog with the primers of the list, one text each
    // (window-modal, returns at once)
    void start(QList<QString> const& primers);

    [[nodiscard]] auto selection() const -> MultiSelection;

private:
    // members are destroyed in reverse order: the widgets before the
    // layouts that hold them
    QVBoxLayout layout_;
    QFormLayout form_;
    std::array<QComboBox, analysis::legacy_multi_primer_fields> fields_;
    QDialogButtonBox buttons_;
};

}  // namespace primerpicker::gui
