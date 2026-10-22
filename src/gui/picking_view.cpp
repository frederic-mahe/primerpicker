#include "gui/picking_view.hpp"

#include "core/analysis/legacy_criteria.hpp"
#include "core/analysis/legacy_pair.hpp"
#include "core/colour/colour.hpp"
#include "core/search/legacy_picking.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"
#include "gui/picking_colours.hpp"

#include <QBrush>
#include <QChar>
#include <QColor>
#include <QEvent>
#include <QFont>
#include <QFontDatabase>
#include <QFontMetrics>
#include <QKeyEvent>
#include <QLabel>
#include <QObject>
#include <QPainter>
#include <QPalette>
#include <QPen>
#include <QPoint>
#include <QPolygon>
#include <QPushButton>
#include <QRect>
#include <QScrollArea>
#include <QSize>
#include <QString>
#include <QWidget>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <optional>
#include <utility>
#include <vector>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// the layout of PRIMER.exe's screen: lines of 60 bases in blocks of 10,
// after the line's first position in 7 columns and 3 spaces
constexpr auto bases_per_line = std::size_t{60};
constexpr auto bases_per_block = std::size_t{10};
constexpr auto number_columns = 7;
constexpr auto gap_columns = 3;
constexpr auto margin = 8;  // pixels around the text
constexpr auto border_width = 2;

auto qcolour(colour::Srgb const colour) -> QColor {
    return QColor{colour.red, colour.green, colour.blue};
}

auto fixed_font() -> QFont {
    return QFontDatabase::systemFont(QFontDatabase::FixedFont);
}

// the first and last positions of the window (top-strand coordinates)
struct Stretch {
    std::size_t first;
    std::size_t last;
};

auto window_stretch(sequence::PrimerLocation const& location) -> Stretch {
    if (location.strand == sequence::Strand::top) {
        return {.first = location.position, .last = location.position + location.length - 1};
    }
    return {.first = location.position - location.length + 1, .last = location.position};
}

auto message_text(search::PickingMessage const message) -> QString {
    switch (message) {
        case search::PickingMessage::primer_too_long:
            return u"Primer Too Long - Extends outside molecule"_s;  // PRIMER.exe's messages
        case search::PickingMessage::primer_must_have_length:
            return u"Primer must have length!"_s;
        case search::PickingMessage::none:
        default:
            return {};
    }
}

// "OK", or the criterion that the dimer reaches
auto dimer_word(bool const rejected, std::size_t const reject_from) -> QString {
    return rejected ? u"! (reject >= %1)"_s.arg(reject_from) : u"OK"_s;
}

}  // namespace

PickingSequence::PickingSequence() {
    setObjectName(u"picking_sequence"_s);
    setFont(fixed_font());
    setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
}

void PickingSequence::show_molecule(std::vector<sequence::Nucleotide> const& molecule,
                                    sequence::Strand const strand) {
    letters_.clear();
    letters_.reserve(molecule.size());
    std::ranges::transform(molecule, std::back_inserter(letters_), [strand](sequence::Nucleotide const nucleotide) {
        return sequence::to_letter(strand == sequence::Strand::top ? nucleotide : sequence::complement(nucleotide));
    });
    updateGeometry();
    resize(sizeHint());
    update();
}

void PickingSequence::show_window(sequence::PrimerLocation const& location, bool const meets_criteria) {
    window_ = location;
    meets_criteria_ = meets_criteria;
    update();
}

auto PickingSequence::sizeHint() const -> QSize {
    auto const metrics = QFontMetrics{font()};
    auto const lines = static_cast<int>((letters_.size() + bases_per_line - 1) / bases_per_line);
    constexpr auto columns =
        number_columns + gap_columns + static_cast<int>(bases_per_line + (bases_per_line / bases_per_block));
    return {(2 * margin) + (columns * metrics.horizontalAdvance(u'M')), (2 * margin) + (lines * 2 * metrics.height())};
}

// a line takes two rows: the 5' mark above, then the bases
auto PickingSequence::base_rectangle(std::size_t const position) const -> QRect {
    assert(position >= 1);
    auto const metrics = QFontMetrics{font()};
    auto const width = metrics.horizontalAdvance(u'M');
    auto const height = metrics.height();
    auto const index = position - 1;
    auto const line = static_cast<int>(index / bases_per_line);
    auto const column = index % bases_per_line;
    auto const text_column = number_columns + gap_columns + static_cast<int>(column + (column / bases_per_block));
    return {margin + (text_column * width), margin + (((2 * line) + 1) * height), width, height};
}

auto PickingSequence::window_rectangle() const -> QRect {
    if (not window_.has_value()) {
        return {};
    }
    auto const [first, last] = window_stretch(*window_);
    return base_rectangle(first).united(base_rectangle(last)).adjusted(0, -base_rectangle(first).height(), 0, 0);
}

auto PickingSequence::palette_colours() const -> PickingPalette const& {
    constexpr auto half = 128;  // of 255
    return palette().color(QPalette::Window).lightness() < half ? dark_picking_palette : light_picking_palette;
}

void PickingSequence::changeEvent(QEvent* const event) {
    if (event->type() == QEvent::PaletteChange or event->type() == QEvent::ApplicationPaletteChange) {
        update();
    }
    QWidget::changeEvent(event);
}

void PickingSequence::paintEvent(QPaintEvent* const /* event */) {
    auto const& colours = palette_colours();
    auto painter = QPainter{this};
    painter.fillRect(rect(), qcolour(colours.view_background));
    painter.setFont(font());
    auto const metrics = QFontMetrics{font()};
    // the window: its background, under the bases
    if (window_.has_value()) {
        auto const [first, last] = window_stretch(*window_);
        auto const fill = qcolour(meets_criteria_ ? colours.meets : colours.fails);
        auto pen = QPen{qcolour(meets_criteria_ ? colours.meets_border : colours.fails_border), border_width,
                        meets_criteria_ ? Qt::SolidLine : Qt::DashLine};
        painter.setPen(pen);
        painter.setBrush(fill);
        for (auto start = first; start <= last;) {
            auto const line_end = std::min(last, ((start - 1) / bases_per_line + 1) * bases_per_line);
            painter.drawRect(base_rectangle(start).united(base_rectangle(line_end)));
            start = line_end + 1;
        }
        // the 5' end: a triangle above it
        auto const five_prime = base_rectangle(window_->position);
        auto const tip = QPoint{five_prime.center().x(), five_prime.top() - 2};
        auto const half_width = five_prime.width() / 2;
        auto const arrow = QPolygon{{tip, QPoint{tip.x() - half_width, tip.y() - half_width},
                                     QPoint{tip.x() + half_width, tip.y() - half_width}}};
        painter.setBrush(pen.color());
        pen.setStyle(Qt::SolidLine);
        painter.setPen(pen);
        painter.drawPolygon(arrow);
    }
    painter.setPen(qcolour(colours.text));
    for (auto start = std::size_t{1}; start <= letters_.size(); start += bases_per_line) {
        auto const first = base_rectangle(start);
        auto const number = QString::number(start).rightJustified(number_columns);
        painter.drawText(QPoint{margin, first.top() + metrics.ascent()}, number);
        auto const end = std::min(letters_.size(), start + bases_per_line - 1);
        for (auto position = start; position <= end; ++position) {
            auto const cell = base_rectangle(position);
            painter.drawText(QPoint{cell.left(), cell.top() + metrics.ascent()},
                             QString{QChar::fromLatin1(letters_.at(position - 1))});
        }
    }
}

PickingView::PickingView()
    : previous_{u"Find &Previous"_s},
      next_{u"Find &Next"_s},
      shorter_{u"Shorter (-)"_s},
      longer_{u"Longer (+)"_s},
      save_{u"&Save to List"_s} {
    setObjectName(u"picking_view"_s);
    setFocusPolicy(Qt::StrongFocus);
    header_.setObjectName(u"picking_header"_s);
    strand_.setObjectName(u"picking_strand"_s);
    verdict_.setObjectName(u"picking_verdict"_s);
    values_.setObjectName(u"picking_values"_s);
    values_.setFont(fixed_font());
    values_.setTextFormat(Qt::PlainText);
    dimers_.setObjectName(u"picking_dimers"_s);
    message_.setObjectName(u"picking_message"_s);
    keys_.setText(u"Move: ← → (one base), ↑ ↓ (one line)    Find next: Ctrl+←, Ctrl+→    "
                  u"Length: + -    Save: S, Enter or Esc"_s);
    keys_.setWordWrap(true);
    scroll_.setObjectName(u"picking_scroll"_s);
    scroll_.setWidget(&sequence_);
    scroll_.setWidgetResizable(false);
    scroll_.setFocusPolicy(Qt::NoFocus);  // the arrows move the primer
    previous_.setObjectName(u"find_previous_button"_s);
    next_.setObjectName(u"find_next_button"_s);
    shorter_.setObjectName(u"shorter_button"_s);
    longer_.setObjectName(u"longer_button"_s);
    save_.setObjectName(u"save_picked_button"_s);
    for (auto* const button : {&previous_, &next_, &shorter_, &longer_, &save_}) {
        button->setFocusPolicy(Qt::NoFocus);
        buttons_layout_.addWidget(button);
    }
    buttons_layout_.addStretch();
    QObject::connect(&previous_, &QPushButton::clicked, this, [this] { find(search::SearchDirection::backward); });
    QObject::connect(&next_, &QPushButton::clicked, this, [this] { find(search::SearchDirection::forward); });
    QObject::connect(&shorter_, &QPushButton::clicked, this, [this] { press(search::PickingKey::shorter); });
    QObject::connect(&longer_, &QPushButton::clicked, this, [this] { press(search::PickingKey::longer); });
    QObject::connect(&save_, &QPushButton::clicked, this, &PickingView::save_requested);
    layout_.addWidget(&header_);
    layout_.addWidget(&strand_);
    layout_.addWidget(&scroll_, 1);
    layout_.addWidget(&verdict_);
    layout_.addWidget(&values_);
    layout_.addWidget(&dimers_);
    layout_.addWidget(&message_);
    layout_.addLayout(&buttons_layout_);
    layout_.addWidget(&keys_);
    setLayout(&layout_);
}

void PickingView::start(Start start) {
    assert(sequence::primer_at(start.molecule, start.location).has_value());
    molecule_name_ = std::move(start.molecule_name);
    molecule_ = std::move(start.molecule);
    location_ = start.location;
    criteria_ = start.criteria;
    partner_ = std::move(start.partner);
    extension_ = std::move(start.extension);
    strand_.setText(location_.strand == sequence::Strand::complement ? u"*** Complementary Strand ***"_s
                                                                    : QString{});
    strand_.setVisible(location_.strand == sequence::Strand::complement);
    dimers_.setVisible(partner_.has_value());
    sequence_.show_molecule(molecule_, location_.strand);
    message_.clear();
    show_state();
    setFocus();
}

auto PickingView::location() const -> sequence::PrimerLocation {
    return location_;
}

void PickingView::press(search::PickingKey const key) {
    auto const [moved, message] = search::legacy_picking_move(location_, molecule_.size(), key);
    show_message(message);
    if (moved.position == location_.position and moved.length == location_.length) {
        return;
    }
    location_ = moved;
    show_state();
    Q_EMIT location_changed();
}

void PickingView::find(search::SearchDirection const direction) {
    auto const found = search::legacy_find_next(molecule_, location_, criteria_, direction, extension_);
    if (not found.has_value()) {
        message_.setText(u"No Primers found"_s);  // PRIMER.exe's message
        return;
    }
    message_.clear();
    location_ = *found;
    show_state();
    Q_EMIT location_changed();
}

void PickingView::show_message(search::PickingMessage const message) {
    message_.setText(message_text(message));
}

void PickingView::show_state() {
    auto const state = search::legacy_picking_state(molecule_, location_, criteria_, extension_);
    assert(state.has_value());
    if (not state.has_value()) {
        return;
    }
    // as PRIMER.exe, the length counts the extension
    // (tests/legacy/extension/picking_extension.tsv)
    auto header = u"Primer length: %1   Molecule: %2"_s.arg(location_.length + extension_.size()).arg(molecule_name_);
    if (not extension_.empty()) {
        auto letters = QString{};
        for (auto const nucleotide : extension_) {
            letters.append(QChar::fromLatin1(sequence::to_letter(nucleotide)));
        }
        header += u"   5'extension: %1"_s.arg(letters);
    }
    header_.setText(header);
    verdict_.setText(state->meets_criteria ? u"Meets criteria selected: Yes"_s : u"Meets criteria selected: No"_s);
    values_.setText(QString::fromStdString(state->value_line));
    sequence_.show_window(location_, state->meets_criteria);
    scroll_.ensureVisible(sequence_.window_rectangle().center().x(), sequence_.window_rectangle().center().y(),
                          sequence_.window_rectangle().width() / 2, sequence_.window_rectangle().height());
    if (partner_.has_value()) {
        auto const primer = sequence::primer_at(molecule_, location_);
        if (primer.has_value()) {
            auto const pair = analysis::PrimerPair{.first = partner_->primer, .second = *primer};
            auto const three_prime = analysis::dimer_at_3prime_end(pair);
            auto const adjacent = analysis::longest_complementary_run(pair);
            auto const flags = analysis::legacy_dimer_flags(three_prime, adjacent, criteria_);
            dimers_.setText(u"A:B dimers with %1: 3' end matches %2 %3, adjacent bases %4 %5"_s.arg(partner_->name)
                                .arg(three_prime)
                                .arg(dimer_word(flags.three_prime, criteria_.reject_dimer_3prime_matches))
                                .arg(adjacent)
                                .arg(dimer_word(flags.adjacent, criteria_.reject_dimer_adjacent_bases)));
        }
    }
}

void PickingView::keyPressEvent(QKeyEvent* const event) {
    auto const control = event->modifiers().testFlag(Qt::ControlModifier);
    switch (event->key()) {
        case Qt::Key_Right:
            if (control) {
                find(search::SearchDirection::forward);
            } else {
                press(search::PickingKey::right);
            }
            return;
        case Qt::Key_Left:
            if (control) {
                find(search::SearchDirection::backward);
            } else {
                press(search::PickingKey::left);
            }
            return;
        case Qt::Key_Down:
            press(search::PickingKey::down);
            return;
        case Qt::Key_Up:
            press(search::PickingKey::up);
            return;
        case Qt::Key_Plus:
            press(search::PickingKey::longer);
            return;
        case Qt::Key_Minus:
            press(search::PickingKey::shorter);
            return;
        case Qt::Key_S:
        case Qt::Key_Return:
        case Qt::Key_Enter:
        case Qt::Key_Escape:
            Q_EMIT save_requested();
            return;
        default:
            QWidget::keyPressEvent(event);
            return;
    }
}

}  // namespace primerpicker::gui
