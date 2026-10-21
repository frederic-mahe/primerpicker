#pragma once

#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractItemModel>
#include <QAction>
#include <QList>
#include <QLineEdit>
#include <QString>
#include <QStringView>
#include <QTableView>
#include <QTest>

namespace primerpicker::gui::test {

// Primer/Create -> Direct Entry of PRIMER, then Enter
inline void direct_entry(MainWindow& window, QString const& primer) {
    child<QAction>(window, u"direct_entry").trigger();
    auto& field = child<QLineEdit>(window, u"direct_entry_field");
    QTest::keyClicks(&field, primer);
    QTest::keyClick(&field, Qt::Key_Return);
}

// Direct Entry of PRIMER, then Save to List
inline void save_direct_entry(MainWindow& window, QString const& primer) {
    direct_entry(window, primer);
    child<QAction>(window, u"save_to_list").trigger();
}

// The texts of the table view named NAME, row by row
inline auto table_rows(MainWindow const& window, QStringView const name) -> QList<QList<QString>> {
    auto const& model = *child<QTableView>(window, name).model();
    auto rows = QList<QList<QString>>{};
    for (auto row = 0; row < model.rowCount(); ++row) {
        auto texts = QList<QString>{};
        for (auto column = 0; column < model.columnCount(); ++column) {
            texts.append(model.data(model.index(row, column)).toString());
        }
        rows.append(texts);
    }
    return rows;
}

}  // namespace primerpicker::gui::test
