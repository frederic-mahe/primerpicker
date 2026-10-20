#pragma once

#include "core/enzymes/enzyme.hpp"
#include "core/enzymes/site_search.hpp"
#include "core/sequence/nucleotide.hpp"
#include "gui/text_table_model.hpp"

#include <QLabel>
#include <QObject>
#include <QTableView>
#include <QVBoxLayout>
#include <QWidget>

#include <span>
#include <vector>

namespace primerpicker::gui {

// The True Sites and Almost Sites screens of PRIMER.exe: the enzyme,
// its recognition sequence and the number of sites, then each site's
// position, strand and (almost sites) match
class SitesView final : public QWidget {
    Q_OBJECT

public:
    SitesView();

    struct Search {
        enzymes::Enzyme const& enzyme;
        enzymes::Mismatches mismatches;
        std::span<sequence::Nucleotide const> molecule;
    };

    void show_sites(Search const& search, std::vector<enzymes::Site> const& sites);

private:
    // members are destroyed in reverse order: the widgets before the
    // layout that holds them
    QVBoxLayout layout_;
    QLabel header_;
    TextTableModel model_;
    QTableView table_;
};

}  // namespace primerpicker::gui
