#include "gui/find_pair_view.hpp"

#include "core/search/legacy_find_pair.hpp"

#include <QAbstractItemView>
#include <QHeaderView>
#include <QList>
#include <QString>

#include <optional>
#include <utility>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// the messages of PRIMER.exe (tests/legacy/README.md, "Limits on
// primers and pairs" and "No Matching Primers Found")
constexpr auto too_many_primers_message = u"More than 100 sites found.  Suggest you redefine matching limits.";
constexpr auto pairs_capped_message = u"More than 1000 possible pairs found!";
constexpr auto no_pair_message = u"No Matching Primers Found";

// a false-priming site Tm, or "-" when there is none
auto site_tm(std::optional<int> const tm_celsius) -> QString {
    return tm_celsius.has_value() ? QString::number(*tm_celsius) : u"-"_s;
}

}  // namespace

FindPairView::FindPairView()
    : model_{{u"Position A"_s, u"Position B"_s, u"GC % A"_s, u"GC % B"_s, u"Tm °C A"_s, u"Tm °C B"_s,
              u"Max 3' dimers"_s, u"Max dimers"_s, u"False priming Tm A"_s, u"False priming Tm B"_s}} {
    setObjectName(u"find_pair_view"_s);
    counts_.setObjectName(u"find_pair_counts"_s);
    message_.setObjectName(u"find_pair_message"_s);
    table_.setObjectName(u"pairs_table"_s);
    table_.setModel(&model_);
    table_.setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_.setSelectionBehavior(QAbstractItemView::SelectRows);
    table_.horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    layout_.addWidget(&counts_);
    layout_.addWidget(&message_);
    layout_.addWidget(&table_);
    setLayout(&layout_);
}

void FindPairView::show_pairs(search::FindPairResult const& result, search::FindPairWeights const& weights) {
    auto rows = QList<QList<QString>>{};
    if (result.limit == search::FindPairLimit::interrupted) {
        counts_.clear();
        message_.setText(u"Search interrupted"_s);
        model_.set_rows(std::move(rows));
        return;
    }
    if (result.limit == search::FindPairLimit::too_many_primers) {
        counts_.clear();
        message_.setText(QString::fromUtf16(too_many_primers_message));
        model_.set_rows(std::move(rows));
        return;
    }
    counts_.setText(u"Primer A: %1 possible primers    Primer B: %2 possible primers    Pairs: %3"_s
                        .arg(result.possible_a)
                        .arg(result.possible_b)
                        .arg(result.pairs.size()));
    if (result.pairs.empty()) {
        message_.setText(QString::fromUtf16(no_pair_message));
    } else if (result.limit == search::FindPairLimit::pairs_capped) {
        message_.setText(QString::fromUtf16(pairs_capped_message));
    } else {
        message_.clear();
    }
    for (auto const& pair : search::legacy_rank_pairs(result.pairs, weights)) {
        rows.append({QString::number(pair.position_a), QString::number(pair.position_b),
                     QString::number(pair.primer_a.gc_percent), QString::number(pair.primer_b.gc_percent),
                     QString::number(pair.primer_a.tm_celsius), QString::number(pair.primer_b.tm_celsius),
                     QString::number(pair.max_dimer_3prime_matches), QString::number(pair.max_dimer_adjacent_bases),
                     site_tm(pair.false_priming_tm_a), site_tm(pair.false_priming_tm_b)});
    }
    model_.set_rows(std::move(rows));
}

}  // namespace primerpicker::gui
