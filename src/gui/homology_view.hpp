#pragma once

#include "core/analysis/legacy_homology.hpp"
#include "core/sequence/nucleotide.hpp"
#include "gui/text_table_model.hpp"

#include <QLabel>
#include <QObject>
#include <QString>
#include <QTableView>
#include <QVBoxLayout>
#include <QWidget>

#include <span>
#include <vector>

namespace primerpicker::gui {

// The Uniqueness and False Priming screens of PRIMER.exe: the sites of
// the molecule where a stretch of the primer matches (its bases, the
// position, the strand and the Tm), or PRIMER.exe's message when there
// is none
class HomologyView final : public QWidget {
    Q_OBJECT

public:
    HomologyView();

    struct Search {
        analysis::HomologyKind kind;
        std::span<sequence::Nucleotide const> primer;
        QString const& molecule_name;
    };

    void show_sites(Search const& search, std::vector<analysis::HomologySite> const& sites);

private:
    // members are destroyed in reverse order: the widgets before the
    // layout that holds them
    QVBoxLayout layout_;
    QLabel header_;
    QLabel message_;
    TextTableModel model_;
    QTableView table_;
};

}  // namespace primerpicker::gui
