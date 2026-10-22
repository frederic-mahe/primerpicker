#include "gui/pair_view.hpp"

#include "core/analysis/legacy_pair.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"

#include <QAbstractItemView>
#include <QHeaderView>
#include <QList>
#include <QString>

#include <cstddef>
#include <span>
#include <utility>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

auto letters(std::span<sequence::Nucleotide const> const nucleotides) -> QString {
    auto text = QString{};
    for (auto const nucleotide : nucleotides) {
        text.append(QChar::fromLatin1(sequence::to_letter(nucleotide)));
    }
    return text;
}

// a value with PRIMER.exe's "!" when flagged
template <typename Value>
auto flagged(Value const value, bool const flag) -> QString {
    return QString::number(value) + (flag ? u"!"_s : QString{});
}

// "N base difference", "N % difference", "N°C difference"; empty when 0
auto difference(int const value, QStringView const unit) -> QString {
    return value == 0 ? QString{} : u"%1%2 difference"_s.arg(value).arg(unit);
}

auto either(bool const first, bool const second, QString const& comment) -> QString {
    return first or second ? comment : QString{};
}

auto hairpin_word(bool const found) -> QString {
    return found ? u"!"_s : u"OK"_s;
}

}  // namespace

PairView::PairView() : model_{{QString{}, u"-A-"_s, u"-B-"_s, u"Comment"_s}} {
    setObjectName(u"pair_view"_s);
    primer_a_.setObjectName(u"pair_a_label"_s);
    primer_b_.setObjectName(u"pair_b_label"_s);
    table_.setObjectName(u"pair_table"_s);
    table_.setModel(&model_);
    table_.setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_.verticalHeader()->hide();
    table_.horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    table_.horizontalHeader()->setStretchLastSection(true);
    layout_.addWidget(&primer_a_);
    layout_.addWidget(&primer_b_);
    layout_.addWidget(&table_);
    setLayout(&layout_);
}

// the whole primers, their 5' extensions included, as PRIMER.exe
// (tests/legacy/extension/ext_pair.screens.txt)
void PairView::show_pair(analysis::LegacyPrimerPair const& primers, analysis::LegacyCriteria const& criteria) {
    primer_a_.setText(u"Primer A: 5'-%1-3'"_s.arg(letters(primers.first.sequence)));
    primer_b_.setText(u"Primer B: 5'-%1-3'"_s.arg(letters(primers.second.sequence)));
    auto const actual_a = analysis::legacy_actual_values(primers.first);
    auto const actual_b = analysis::legacy_actual_values(primers.second);
    auto const cautions_a = analysis::legacy_pair_cautions(actual_a, criteria);
    auto const cautions_b = analysis::legacy_pair_cautions(actual_b, criteria);
    auto const pair = analysis::legacy_pair_values(primers);
    auto const flags = [&criteria](std::size_t const three_prime, std::size_t const adjacent) {
        return analysis::legacy_dimer_flags(three_prime, adjacent, criteria);
    };
    auto const flags_aa = flags(actual_a.dimer_3prime_matches, actual_a.dimer_adjacent_bases);
    auto const flags_bb = flags(actual_b.dimer_3prime_matches, actual_b.dimer_adjacent_bases);
    auto const flags_ab = flags(pair.dimer_3prime_matches, pair.dimer_adjacent_bases);
    auto rows = QList<QList<QString>>{
        {u"Length"_s, QString::number(primers.first.sequence.size()), QString::number(primers.second.sequence.size()),
         difference(pair.length_difference, u" base")},
        {u"% GC"_s, flagged(actual_a.gc_percent, cautions_a.gc), flagged(actual_b.gc_percent, cautions_b.gc),
         difference(pair.gc_difference_percent, u" %")},
        {u"Tm °C"_s, flagged(actual_a.tm_celsius, cautions_a.tm), flagged(actual_b.tm_celsius, cautions_b.tm),
         difference(pair.tm_difference_celsius, u"°C")},
        {u"Runs of G/C at 3' end"_s, flagged(actual_a.gc_run_3prime, cautions_a.gc_run_3prime),
         flagged(actual_b.gc_run_3prime, cautions_b.gc_run_3prime),
         either(cautions_a.gc_run_3prime, cautions_b.gc_run_3prime, u"Runs found"_s)},
        {u"Runs of bases (any)"_s, flagged(actual_a.longest_run, cautions_a.runs),
         flagged(actual_b.longest_run, cautions_b.runs), either(cautions_a.runs, cautions_b.runs, u"Runs found"_s)},
        {u"Hairpins"_s, hairpin_word(cautions_a.hairpin), hairpin_word(cautions_b.hairpin),
         either(cautions_a.hairpin, cautions_b.hairpin, u"Hairpin(s) found"_s)},
        // the dimers: A:A, B:B and A:B in the three columns
        {u"3' end matches"_s, flagged(actual_a.dimer_3prime_matches, flags_aa.three_prime),
         flagged(actual_b.dimer_3prime_matches, flags_bb.three_prime),
         flagged(pair.dimer_3prime_matches, flags_ab.three_prime)},
        {u"Any adj homol bases"_s, flagged(actual_a.dimer_adjacent_bases, flags_aa.adjacent),
         flagged(actual_b.dimer_adjacent_bases, flags_bb.adjacent),
         flagged(pair.dimer_adjacent_bases, flags_ab.adjacent)},
    };
    model_.set_rows(std::move(rows));
}

}  // namespace primerpicker::gui
