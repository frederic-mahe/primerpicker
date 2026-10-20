#include "gui/sites_view.hpp"

#include "core/enzymes/enzyme.hpp"
#include "core/enzymes/site_search.hpp"
#include "core/sequence/iupac.hpp"
#include "core/sequence/strand.hpp"

#include <QAbstractItemView>
#include <QHeaderView>
#include <QList>
#include <QString>

#include <utility>
#include <vector>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

SitesView::SitesView() : model_{{u"Position"_s, u"Strand"_s, u"Match"_s}} {
    setObjectName(u"sites_view"_s);
    header_.setObjectName(u"sites_header"_s);
    table_.setObjectName(u"sites_table"_s);
    table_.setModel(&model_);
    table_.setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_.verticalHeader()->hide();
    table_.horizontalHeader()->setStretchLastSection(true);
    layout_.addWidget(&header_);
    layout_.addWidget(&table_);
    setLayout(&layout_);
}

void SitesView::show_sites(Search const& search, std::vector<enzymes::Site> const& sites) {
    auto const name = QString::fromLatin1(search.enzyme.name);
    auto const site = QString::fromLatin1(search.enzyme.site);
    auto const almost = search.mismatches == enzymes::Mismatches::at_most_one;
    // as PRIMER.exe, the Almost Sites header shows the reverse complement
    // of a non-palindromic site
    auto const reversed = sequence::reverse_complement_iupac(search.enzyme.site);
    auto const complement = almost and reversed.has_value() and *reversed != search.enzyme.site
                                ? u" (%1)"_s.arg(QString::fromStdString(*reversed))
                                : QString{};
    header_.setText(u"Search for: %1 %2%3    Number of Sites = %4"_s.arg(name, site, complement)
                        .arg(sites.size()));
    auto rows = QList<QList<QString>>{};
    for (auto const& found : sites) {
        auto const strand = found.strand == sequence::Strand::complement ? u"complement"_s : u"top"_s;
        auto const match =
            almost ? QString::fromStdString(enzymes::site_match_text(found, search.molecule, search.enzyme.site.size()))
                   : QString{};
        rows.append({QString::number(found.position), strand, match});
    }
    model_.set_rows(std::move(rows));
}

}  // namespace primerpicker::gui
