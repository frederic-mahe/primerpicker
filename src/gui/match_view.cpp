#include "gui/match_view.hpp"

#include "core/search/legacy_match.hpp"

#include <QAbstractItemView>
#include <QHeaderView>
#include <QList>
#include <QString>

#include <utility>
#include <vector>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// PRIMER.exe's box when no primer B matches (tests/legacy/README.md,
// "The Match Active Primer report")
constexpr auto no_match_message = u"No Matching Primers Found";

}  // namespace

MatchView::MatchView()
    : model_{{u"Position"_s, u"bps from target"_s, u"GC %"_s, u"Tm °C"_s, u"Self 3' dimers"_s, u"Self dimers"_s,
              u"A:B 3' dimers"_s, u"A:B dimers"_s, u"Runs"_s, u"3' G/C runs"_s}} {
    setObjectName(u"match_view"_s);
    header_.setObjectName(u"match_header"_s);
    message_.setObjectName(u"match_message"_s);
    table_.setObjectName(u"match_table"_s);
    table_.setModel(&model_);
    table_.setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_.setSelectionBehavior(QAbstractItemView::SelectRows);
    table_.verticalHeader()->hide();
    table_.horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    layout_.addWidget(&header_);
    layout_.addWidget(&message_);
    layout_.addWidget(&table_);
    setLayout(&layout_);
}

void MatchView::show_list(QString const& header, std::vector<search::MatchCandidate> const& list) {
    header_.setText(header);
    if (list.empty()) {
        message_.setText(QString::fromUtf16(no_match_message));
    } else {
        message_.clear();
    }
    auto rows = QList<QList<QString>>{};
    for (auto const& candidate : list) {
        auto const& actual = candidate.actual;
        rows.append({QString::number(candidate.position), QString::number(candidate.distance),
                     QString::number(actual.gc_percent), QString::number(actual.tm_celsius),
                     QString::number(actual.dimer_3prime_matches), QString::number(actual.dimer_adjacent_bases),
                     QString::number(candidate.inter_dimer_3prime_matches),
                     QString::number(candidate.inter_dimer_adjacent_bases), QString::number(actual.longest_run),
                     QString::number(actual.gc_run_3prime)});
    }
    model_.set_rows(std::move(rows));
}

}  // namespace primerpicker::gui
