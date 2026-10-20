#include "gui/analysis_view.hpp"

#include "core/analysis/legacy_hairpin.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/composition.hpp"
#include "core/sequence/nucleotide.hpp"

#include <QAbstractItemView>
#include <QAbstractScrollArea>
#include <QFontDatabase>
#include <QHeaderView>
#include <QLabel>
#include <QList>
#include <QObject>
#include <QPushButton>
#include <QSizePolicy>
#include <QString>

#include <algorithm>
#include <array>
#include <cstddef>
#include <format>
#include <iterator>
#include <span>
#include <string>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

auto from_utf8(std::string const& text) -> QString {
    return QString::fromStdString(text);
}

auto verdict_word(bool const meets) -> QString {
    return meets ? u"YES"_s : u"NO"_s;
}

auto letters(std::span<sequence::Nucleotide const> const nucleotides) -> std::string {
    auto text = std::string{};
    std::ranges::transform(nucleotides, std::back_inserter(text), sequence::to_letter);
    return text;
}

// A row of the criteria table: criterion, setting, actual value and
// whether the primer meets it
struct CriterionRow {
    QString criterion;
    QString setting;
    QString actual;
    bool meets = false;
};

// The labels and settings of PRIMER.exe's One Primer report
// (tests/legacy/reports/INSPECT.PRN)
auto criteria_rows(analysis::LegacyActualValues const& actual, analysis::LegacyCriteria const& criteria,
                   analysis::LegacyVerdict const& verdict) -> QList<QList<QString>> {
    auto const count = [](std::size_t const value) { return QString::number(value); };
    constexpr auto criterion_count = std::size_t{7};
    auto const rows = std::array<CriterionRow, criterion_count>{{
        {.criterion = u"% GC"_s,
         .setting = from_utf8(std::format("Min {}, Max {}", criteria.gc_min_percent, criteria.gc_max_percent)),
         .actual = QString::number(actual.gc_percent),
         .meets = verdict.gc},
        {.criterion = u"Tm °C"_s,
         .setting = from_utf8(std::format("Min {}, Max {}", criteria.tm_min_celsius, criteria.tm_max_celsius)),
         .actual = QString::number(actual.tm_celsius),
         .meets = verdict.tm},
        // as on the Primer Analysis screen, the hairpin energy is not in
        // this column (see the hairpin line)
        {.criterion = u"No Hairpins"_s,
         .setting = from_utf8(std::format("Energy cutoff {:.1f} kcal", criteria.hairpin_cutoff.kcal_per_mol)),
         .actual = u"--"_s,
         .meets = verdict.no_hairpin},
        {.criterion = u"No 3' Dimers"_s,
         .setting = from_utf8(std::format("Reject >= {} matches at 3' end", criteria.reject_dimer_3prime_matches)),
         .actual = count(actual.dimer_3prime_matches),
         .meets = verdict.no_3prime_dimer},
        {.criterion = u"No Dimers"_s,
         .setting = from_utf8(std::format("Reject >= {} adj homol bases", criteria.reject_dimer_adjacent_bases)),
         .actual = count(actual.dimer_adjacent_bases),
         .meets = verdict.no_dimer},
        {.criterion = u"No Runs"_s,
         .setting = from_utf8(std::format("Reject >= {} base runs", criteria.reject_run_length)),
         .actual = count(actual.longest_run),
         .meets = verdict.no_runs},
        {.criterion = u"No 3'GC runs"_s,
         .setting = from_utf8(std::format("Reject >= {} G or C at 3' end", criteria.reject_gc_run_3prime)),
         .actual = count(actual.gc_run_3prime),
         .meets = verdict.no_3prime_gc_run},
    }};
    auto table = QList<QList<QString>>{};
    for (auto const& row : rows) {
        table.append({row.criterion, row.setting, row.actual, verdict_word(row.meets)});
    }
    return table;
}

auto hairpin_text(analysis::LegacyActualValues const& actual) -> QString {
    if (not actual.hairpin.has_value()) {
        return u"No Hairpins Found"_s;
    }
    return from_utf8(std::format("{:.1f} kcal", actual.hairpin->kcal_per_mol));
}

// The drawing of the Hairpins screen (H), with the screen's columns, or
// its "No Hairpins Found"; the whole sequence, 5' extension included
auto hairpin_drawing_text(std::span<sequence::Nucleotide const> const sequence) -> QString {
    auto const hairpin = analysis::legacy_hairpin(sequence);
    if (not hairpin.has_value()) {
        return u"No Hairpins Found"_s;
    }
    auto const drawing = analysis::legacy_hairpin_drawing(sequence, *hairpin);
    return from_utf8(std::format("{}\n{}\n{}", drawing.top, drawing.marks, drawing.bottom));
}

}  // namespace

AnalysisView::AnalysisView()
    : criteria_model_{{u"Criteria"_s, u"Criteria Setting"_s, u"Actual"_s, u"Meets Criteria"_s}},
      print_{u"&Print to Disk..."_s},
      uniqueness_{u"&Uniqueness"_s},
      false_priming_{u"&False Priming"_s} {
    setObjectName(u"analysis_view"_s);
    sequence_.setObjectName(u"sequence_label"_s);
    length_.setObjectName(u"length_label"_s);
    extension_.setObjectName(u"extension_label"_s);
    composition_.setObjectName(u"composition_label"_s);
    hairpin_.setObjectName(u"hairpin_label"_s);
    hairpin_drawing_.setObjectName(u"hairpin_drawing"_s);
    hairpin_drawing_.setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    hairpin_drawing_.setTextFormat(Qt::PlainText);
    verdict_.setObjectName(u"verdict_label"_s);
    criteria_table_.setObjectName(u"criteria_table"_s);
    sequence_.setTextInteractionFlags(Qt::TextSelectableByMouse);
    criteria_table_.setModel(&criteria_model_);
    criteria_table_.setEditTriggers(QAbstractItemView::NoEditTriggers);
    criteria_table_.setSelectionMode(QAbstractItemView::NoSelection);
    criteria_table_.verticalHeader()->hide();
    criteria_table_.horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    criteria_table_.horizontalHeader()->setStretchLastSection(true);
    // as tall as its rows
    criteria_table_.setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContents);
    criteria_table_.setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    criteria_table_.setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Maximum);

    primer_layout_.addRow(u"Sequence:"_s, &sequence_);
    primer_layout_.addRow(u"Length:"_s, &length_);
    primer_layout_.addRow(u"5' Extension:"_s, &extension_);
    primer_layout_.addRow(u"Composition:"_s, &composition_);
    result_layout_.addRow(u"Hairpin:"_s, &hairpin_);
    result_layout_.addRow(&hairpin_drawing_);
    result_layout_.addRow(u"Meets criteria selected:"_s, &verdict_);
    layout_.addLayout(&primer_layout_);
    layout_.addWidget(&criteria_table_);
    layout_.addLayout(&result_layout_);
    print_.setObjectName(u"print_button"_s);
    uniqueness_.setObjectName(u"uniqueness_button"_s);
    false_priming_.setObjectName(u"false_priming_button"_s);
    QObject::connect(&print_, &QPushButton::clicked, this, &AnalysisView::print_requested);
    QObject::connect(&uniqueness_, &QPushButton::clicked, this, &AnalysisView::uniqueness_requested);
    QObject::connect(&false_priming_, &QPushButton::clicked, this, &AnalysisView::false_priming_requested);
    buttons_layout_.addWidget(&print_);
    buttons_layout_.addWidget(&uniqueness_);
    buttons_layout_.addWidget(&false_priming_);
    buttons_layout_.addStretch();
    layout_.addLayout(&buttons_layout_);
    layout_.addStretch();
    setLayout(&layout_);
}

void AnalysisView::show_primer(analysis::LegacyPrimer const& primer, analysis::LegacyCriteria const& criteria) {
    auto const actual = analysis::legacy_actual_values(primer);
    auto const verdict = analysis::evaluate(actual, criteria);
    auto const extension = primer.sequence.first(primer.extension_length);
    auto const counts = sequence::composition(primer.sequence);
    // as the One Primer report: 5'-EXTENSION-PRIMER-3'
    auto const without_extension = primer.sequence.subspan(primer.extension_length);
    sequence_.setText(from_utf8(extension.empty()
                                    ? std::format("5'-{}-3'", letters(without_extension))
                                    : std::format("5'-{}-{}-3'", letters(extension), letters(without_extension))));
    length_.setText(QString::number(primer.sequence.size()));
    extension_.setText(extension.empty() ? u"--none--"_s : from_utf8(letters(extension)));
    composition_.setText(from_utf8(std::format("A: {}  C: {}  G: {}  T: {}", counts.adenine, counts.cytosine,
                                               counts.guanine, counts.thymine)));
    criteria_model_.set_rows(criteria_rows(actual, criteria, verdict));
    hairpin_.setText(hairpin_text(actual));
    hairpin_drawing_.setText(hairpin_drawing_text(primer.sequence));
    verdict_.setText(verdict_word(verdict.all()));
}

}  // namespace primerpicker::gui
