#pragma once

#include "core/search/legacy_find_pair.hpp"
#include "gui/text_table_model.hpp"

#include <QLabel>
#include <QObject>
#include <QTableView>
#include <QVBoxLayout>
#include <QWidget>

namespace primerpicker::gui {

// The Find Primer Pair list of PRIMER.exe: the numbers of possible
// primers and pairs, PRIMER.exe's message when the search stops or
// finds nothing, and the pairs ranked with the weighting factors
class FindPairView final : public QWidget {
    Q_OBJECT

public:
    FindPairView();

    void show_pairs(search::FindPairResult const& result, search::FindPairWeights const& weights);

private:
    // members are destroyed in reverse order: the widgets before the
    // layout that holds them
    QVBoxLayout layout_;
    QLabel counts_;
    QLabel message_;
    TextTableModel model_;
    QTableView table_;
};

}  // namespace primerpicker::gui
