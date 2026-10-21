#pragma once

#include "core/search/legacy_match.hpp"
#include "gui/text_table_model.hpp"

#include <QLabel>
#include <QObject>
#include <QString>
#include <QTableView>
#include <QVBoxLayout>
#include <QWidget>

#include <vector>

namespace primerpicker::gui {

// The Find Matching Primer list of PRIMER.exe: the primers B that match
// the active primer, by position (distance to the target, GC %, Tm,
// self-dimers, A:B dimers, runs, 3' G/C runs), or PRIMER.exe's "No
// Matching Primers Found"
class MatchView final : public QWidget {
    Q_OBJECT

public:
    MatchView();

    void show_list(QString const& header, std::vector<search::MatchCandidate> const& list);

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
