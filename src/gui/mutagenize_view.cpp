#include "gui/mutagenize_view.hpp"

#include "core/analysis/legacy_mutagenize.hpp"
#include "core/enzymes/rebase_enzymes.hpp"
#include "core/search/legacy_picking.hpp"
#include "core/sequence/codon_table.hpp"
#include "core/sequence/genetic_code.hpp"
#include "core/sequence/nucleotide.hpp"

#include <QAbstractItemView>
#include <QDialog>
#include <QEvent>
#include <QFontDatabase>
#include <QKeyEvent>
#include <QLabel>
#include <QList>
#include <QListWidget>
#include <QObject>
#include <QPushButton>
#include <QString>
#include <QWidget>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

constexpr auto codon_length = std::size_t{3};

auto from_utf8(std::string const& text) -> QString {
    return QString::fromStdString(text);
}

auto base_of_key(int const key) -> std::optional<sequence::Nucleotide> {
    switch (key) {
        case Qt::Key_A:
            return sequence::Nucleotide::adenine;
        case Qt::Key_C:
            return sequence::Nucleotide::cytosine;
        case Qt::Key_G:
            return sequence::Nucleotide::guanine;
        case Qt::Key_T:
            return sequence::Nucleotide::thymine;
        default:
            return std::nullopt;
    }
}

// The lines of the screen as plain text, the base under the edit
// highlight in bold and reverse video (HTML)
auto screen_html(analysis::MutagenizeView const& view, std::size_t const cursor) -> QString {
    auto const escaped = [](std::string const& line) { return from_utf8(line).toHtmlEscaped(); };
    auto lines = QList<QString>{};
    for (auto const& frame : view.original_frames) {
        lines.append(escaped(frame));
    }
    lines.append(escaped(view.molecule));
    auto const column = analysis::detail::mutagenize_primer_column + cursor;
    auto const& primer = view.primer;
    auto const before = primer.substr(0, std::min(column, primer.size()));
    auto const under = column < primer.size() ? primer.substr(column, 1) : std::string{" "};
    auto const after = column + 1 < primer.size() ? primer.substr(column + 1) : std::string{};
    lines.append(escaped(before) + u"<span style=\"background-color: palette(highlight); "
                                   u"color: palette(highlighted-text); font-weight: bold\">"_s +
                 escaped(under) + u"</span>"_s + escaped(after));
    for (auto const& frame : view.mutated_frames) {
        lines.append(escaped(frame));
    }
    return u"<pre>"_s + lines.join(u'\n') + u"</pre>"_s;
}

}  // namespace

MutagenizeListDialog::MutagenizeListDialog(QWidget& parent, ListNames const& names)
    : QDialog{&parent}, help_{u"Type an initial letter to jump to it    Esc = Exit"_s} {
    setObjectName(names.object_name);
    setWindowTitle(names.title);
    list_.setObjectName(names.object_name + u"_list"_s);
    list_.setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    list_.setSelectionMode(QAbstractItemView::SingleSelection);
    layout_.addWidget(&list_);
    layout_.addWidget(&help_);
    setLayout(&layout_);
}

void MutagenizeListDialog::start(QList<QString> const& rows, int const row) {
    list_.clear();
    list_.addItems(rows);
    list_.setCurrentRow(row);
    list_.setFocus();
    open();
}

MutagenizeView::MutagenizeView()
    : enzymes_button_{u"&Enzymes (F1)"_s},
      codons_button_{u"&Codons (F2)"_s},
      inspect_button_{u"&Inspect"_s},
      save_button_{u"&Save"_s},
      enzymes_{*this, {.title = u"Enzymes"_s, .object_name = u"mutagenize_enzymes"_s}},
      codons_{*this, {.title = u"Codon Table"_s, .object_name = u"mutagenize_codons"_s}} {
    setObjectName(u"mutagenize_view"_s);
    setFocusPolicy(Qt::StrongFocus);
    header_.setObjectName(u"mutagenize_header"_s);
    screen_.setObjectName(u"mutagenize_screen"_s);
    screen_.setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    screen_.setTextFormat(Qt::RichText);
    original_values_.setObjectName(u"mutagenize_original"_s);
    mutated_values_.setObjectName(u"mutagenize_mutated"_s);
    for (auto* const label : {&original_values_, &mutated_values_}) {
        label->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
        label->setTextFormat(Qt::PlainText);
    }
    message_.setObjectName(u"mutagenize_message"_s);
    keys_.setText(u"Edit: A C G T over the highlighted base    Move: ← →    F1: enzymes    F2: codons    "
                  u"I: inspect    S: save    Esc: exit"_s);
    keys_.setWordWrap(true);
    enzymes_button_.setObjectName(u"mutagenize_enzymes_button"_s);
    codons_button_.setObjectName(u"mutagenize_codons_button"_s);
    inspect_button_.setObjectName(u"mutagenize_inspect_button"_s);
    save_button_.setObjectName(u"mutagenize_save_button"_s);
    for (auto* const button : {&enzymes_button_, &codons_button_, &inspect_button_, &save_button_}) {
        button->setFocusPolicy(Qt::NoFocus);
        buttons_layout_.addWidget(button);
    }
    buttons_layout_.addStretch();
    QObject::connect(&enzymes_button_, &QPushButton::clicked, this, &MutagenizeView::show_enzymes);
    QObject::connect(&codons_button_, &QPushButton::clicked, this, &MutagenizeView::show_codons);
    QObject::connect(&inspect_button_, &QPushButton::clicked, this, &MutagenizeView::inspect_requested);
    QObject::connect(&save_button_, &QPushButton::clicked, this, &MutagenizeView::save_requested);
    layout_.addWidget(&header_);
    layout_.addWidget(&screen_);
    layout_.addWidget(&original_values_);
    layout_.addWidget(&mutated_values_);
    layout_.addWidget(&message_);
    layout_.addLayout(&buttons_layout_);
    layout_.addWidget(&keys_);
    layout_.addStretch();
    setLayout(&layout_);
}

void MutagenizeView::start(Start start) {
    assert(start.position >= 1 and start.position - 1 + start.primer.size() <= start.molecule.size());
    header_.setText(u"Mutagenize Primer    Molecule: %1    Position: %2"_s.arg(start.molecule_name)
                        .arg(start.position));
    molecule_ = std::move(start.molecule);
    position_ = start.position;
    original_ = std::move(start.primer);
    mutated_ = original_;
    cursor_ = 0;
    criteria_ = start.criteria;
    message_.clear();
    show_state();
    setFocus();
}

auto MutagenizeView::original() const -> std::vector<sequence::Nucleotide> const& {
    return original_;
}

auto MutagenizeView::mutated() const -> std::vector<sequence::Nucleotide> const& {
    return mutated_;
}

auto MutagenizeView::cursor() const -> std::size_t {
    return cursor_;
}

void MutagenizeView::move_cursor(int const bases) {
    if (mutated_.empty()) {
        return;
    }
    auto const last = static_cast<int>(mutated_.size()) - 1;
    cursor_ = static_cast<std::size_t>(std::clamp(static_cast<int>(cursor_) + bases, 0, last));
    message_.clear();
    show_state();
}

void MutagenizeView::type_base(sequence::Nucleotide const base) {
    if (cursor_ >= mutated_.size()) {
        return;
    }
    mutated_.at(cursor_) = base;  // overwrite; the highlight stays
    message_.clear();
    show_state();
}

void MutagenizeView::show_saved() {
    message_.setText(u"Primer Saved"_s);  // PRIMER.exe's message
}

void MutagenizeView::show_state() {
    auto const view = analysis::legacy_mutagenize_view(molecule_, position_, mutated_);
    screen_.setText(screen_html(view, cursor_));
    original_values_.setText(u"ORIGINAL: "_s + from_utf8(search::legacy_mutagenize_line(original_, criteria_)));
    mutated_values_.setText(u"MUTATED:  "_s + from_utf8(search::legacy_mutagenize_line(mutated_, criteria_)));
}

// the enzymes PrimerPicker ships, with their sites (PRIMER.exe lists its
// own enzymes; choosing one changes nothing)
void MutagenizeView::show_enzymes() {
    auto rows = QList<QString>{};
    for (auto const& enzyme : enzymes::rebase_enzymes) {
        constexpr auto name_width = -12;  // left-aligned in 12 columns
        rows.append(u"%1 %2"_s.arg(QString::fromLatin1(enzyme.name), name_width).arg(QString::fromLatin1(enzyme.site)));
    }
    enzymes_.start(rows, 0);
}

// opened on the amino acid of the codon (first frame of the mutated
// primer, M) that holds the highlighted base
void MutagenizeView::show_codons() {
    auto rows = QList<QString>{};
    auto current = 0;
    auto const codon_start = cursor_ - (cursor_ % codon_length);
    auto const code = codon_start + codon_length <= mutated_.size()
                          ? sequence::translate(mutated_.at(codon_start), mutated_.at(codon_start + 1),
                                                mutated_.at(codon_start + 2))
                          : ' ';
    for (auto const& entry : sequence::legacy_codon_table) {
        if (entry.code == code) {
            current = static_cast<int>(rows.size());
        }
        rows.append(u"%1  %2  %3"_s.arg(QString::fromLatin1(entry.name))
                        .arg(QChar::fromLatin1(entry.code))
                        .arg(QString::fromLatin1(entry.codons)));
    }
    codons_.start(rows, current);
}

auto MutagenizeView::event(QEvent* const event) -> bool {
    auto const* const key_event = event->type() == QEvent::ShortcutOverride ? dynamic_cast<QKeyEvent*>(event) : nullptr;
    if (key_event != nullptr and (key_event->key() == Qt::Key_F1 or key_event->key() == Qt::Key_F2)) {
        event->accept();
        return true;
    }
    return QWidget::event(event);
}

void MutagenizeView::keyPressEvent(QKeyEvent* const event) {
    if (auto const base = base_of_key(event->key()); base.has_value()) {
        type_base(*base);
        return;
    }
    switch (event->key()) {
        case Qt::Key_Right:
            move_cursor(1);
            return;
        case Qt::Key_Left:
            move_cursor(-1);
            return;
        case Qt::Key_F1:
            show_enzymes();
            return;
        case Qt::Key_F2:
            show_codons();
            return;
        case Qt::Key_I:
            Q_EMIT inspect_requested();
            return;
        case Qt::Key_S:
            Q_EMIT save_requested();
            return;
        case Qt::Key_Escape:
            Q_EMIT exit_requested();
            return;
        default:
            QWidget::keyPressEvent(event);
            return;
    }
}

}  // namespace primerpicker::gui
