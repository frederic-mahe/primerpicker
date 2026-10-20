#include "gui/homology_view.hpp"

#include "core/analysis/legacy_homology.hpp"
#include "core/sequence/nucleotide.hpp"
#include "core/sequence/strand.hpp"

#include <QAbstractItemView>
#include <QHeaderView>
#include <QList>
#include <QString>

#include <span>
#include <utility>
#include <vector>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// the messages of PRIMER.exe when there is no site (tests/legacy/README.md,
// "Uniqueness and false priming")
constexpr auto unique_message = u"Primer is unique";
constexpr auto no_false_priming_message = u"No false priming found";

auto letters(std::span<sequence::Nucleotide const> const nucleotides) -> QString {
    auto text = QString{};
    for (auto const nucleotide : nucleotides) {
        text.append(QChar::fromLatin1(sequence::to_letter(nucleotide)));
    }
    return text;
}

}  // namespace

HomologyView::HomologyView() : model_{{u"Homologous bases"_s, u"Position"_s, u"Strand"_s, u"Tm °C"_s}} {
    setObjectName(u"homology_view"_s);
    header_.setObjectName(u"homology_header"_s);
    message_.setObjectName(u"homology_message"_s);
    table_.setObjectName(u"homology_table"_s);
    table_.setModel(&model_);
    table_.setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_.verticalHeader()->hide();
    table_.horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    layout_.addWidget(&header_);
    layout_.addWidget(&message_);
    layout_.addWidget(&table_);
    setLayout(&layout_);
}

void HomologyView::show_sites(Search const& search, std::vector<analysis::HomologySite> const& sites) {
    auto const uniqueness = search.kind == analysis::HomologyKind::uniqueness;
    header_.setText(u"%1 of 5'-%2-3' on %3"_s.arg(uniqueness ? u"Uniqueness"_s : u"False Priming"_s,
                                                   letters(search.primer), search.molecule_name));
    if (sites.empty()) {
        message_.setText(QString::fromUtf16(uniqueness ? unique_message : no_false_priming_message));
    } else {
        message_.clear();
    }
    auto rows = QList<QList<QString>>{};
    for (auto const& site : sites) {
        auto const strand = site.strand == sequence::Strand::complement ? u"complement"_s : u"top"_s;
        rows.append({letters(search.primer.subspan(site.primer_start, site.length)), QString::number(site.position),
                     strand, QString::number(site.tm_celsius)});
    }
    model_.set_rows(std::move(rows));
}

}  // namespace primerpicker::gui
