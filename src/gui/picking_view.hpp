#pragma once

#include "core/analysis/legacy_verdict.hpp"
#include "core/search/legacy_picking.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "gui/picking_colours.hpp"

#include <QEvent>
#include <QHBoxLayout>
#include <QKeyEvent>
#include <QLabel>
#include <QObject>
#include <QPaintEvent>
#include <QPushButton>
#include <QRect>
#include <QScrollArea>
#include <QSize>
#include <QString>
#include <QVBoxLayout>
#include <QWidget>

#include <cstddef>
#include <optional>
#include <vector>

namespace primerpicker::gui {

// The molecule as the picking screen shows it: lines of 60 bases in
// blocks of 10, each line after its first position, the bases of the
// primer's strand (complemented for the complement strand, as
// PRIMER.exe), and the primer's window: its background and border tell
// whether it meets the criteria (picking_colours.hpp), a mark above
// the line shows its 5' end
class PickingSequence final : public QWidget {
public:
    PickingSequence();

    void show_molecule(std::vector<sequence::Nucleotide> const& molecule, sequence::Strand strand);
    void show_window(sequence::PrimerLocation const& location, bool meets_criteria);
    // where the window is drawn, to scroll to it
    [[nodiscard]] auto window_rectangle() const -> QRect;
    [[nodiscard]] auto sizeHint() const -> QSize override;

protected:
    void paintEvent(QPaintEvent* event) override;
    void changeEvent(QEvent* event) override;

private:
    [[nodiscard]] auto palette_colours() const -> PickingPalette const&;
    [[nodiscard]] auto base_rectangle(std::size_t position) const -> QRect;

    std::vector<char> letters_;
    std::optional<sequence::PrimerLocation> window_;
    bool meets_criteria_ = false;
};

// The live state of the primer picked second, for the A:B dimers with
// the primer picked first (not in PRIMER.exe)
struct PickingPartner {
    std::vector<sequence::Nucleotide> primer;  // without 5' extension
    QString name;                               // as the Primer List shows it
};

// Primer/Create -> From Molecule and Modify -> Position: PRIMER.exe's
// MODIFY PRIMER POSITION screen. The arrows move the primer by one
// base, Up and Down by a line of 60 bases, + and - change its length,
// Ctrl+Right and Ctrl+Left find the next or previous primer that meets
// the criteria; S, Enter or Esc save it to the Primer List (Esc too, as
// in PRIMER.exe). The values are updated on every move. With a partner,
// a separate line shows the A:B dimers with it.
class PickingView final : public QWidget {
    Q_OBJECT

public:
    PickingView();

    struct Start {
        QString molecule_name;
        std::vector<sequence::Nucleotide> molecule;
        sequence::PrimerLocation location;  // of the molecule's bases
        analysis::LegacyCriteria criteria;
        std::optional<PickingPartner> partner;
        std::vector<sequence::Nucleotide> extension = {};  // the primer's 5' extension
    };
    // Requires a location that fits on the molecule
    void start(Start start);

    [[nodiscard]] auto location() const -> sequence::PrimerLocation;

    // the keys of the screen, also behind its buttons
    void press(search::PickingKey key);
    void find(search::SearchDirection direction);

Q_SIGNALS:
    // the primer moved or changed length
    void location_changed();
    // S, Enter or Esc
    void save_requested();

protected:
    void keyPressEvent(QKeyEvent* event) override;

private:
    void show_state();
    void show_message(search::PickingMessage message);

    QString molecule_name_;
    std::vector<sequence::Nucleotide> molecule_;
    sequence::PrimerLocation location_{.position = 1, .length = 1, .strand = sequence::Strand::top};
    analysis::LegacyCriteria criteria_;
    std::vector<sequence::Nucleotide> extension_;
    std::optional<PickingPartner> partner_;
    // members are destroyed in reverse order: the widgets before the
    // layouts that hold them
    QVBoxLayout layout_;
    QHBoxLayout buttons_layout_;
    QLabel header_;
    QLabel strand_;
    QScrollArea scroll_;
    // after the scroll area that shows it: destroyed first, it leaves
    // the scroll area, which would otherwise delete it
    PickingSequence sequence_;
    QLabel verdict_;
    QLabel values_;
    QLabel dimers_;
    QLabel message_;
    QLabel keys_;
    QPushButton previous_;
    QPushButton next_;
    QPushButton shorter_;
    QPushButton longer_;
    QPushButton save_;
};

}  // namespace primerpicker::gui
