#pragma once

#include <QAbstractTableModel>
#include <QList>
#include <QModelIndex>
#include <QObject>
#include <QString>
#include <QVariant>

namespace primerpicker::gui {

// A read-only table of texts (the criteria of a primer, the sites of an
// enzyme), as a Qt model: column headers, then rows of as many texts
class TextTableModel final : public QAbstractTableModel {
    Q_OBJECT

public:
    explicit TextTableModel(QList<QString> headers);

    void set_rows(QList<QList<QString>> rows);

    // a tooltip per row (the same for every cell of the row)
    void set_tooltips(QList<QString> tooltips);

    [[nodiscard]] auto rowCount(QModelIndex const& parent) const -> int override;
    [[nodiscard]] auto columnCount(QModelIndex const& parent) const -> int override;
    [[nodiscard]] auto data(QModelIndex const& index, int role) const -> QVariant override;
    [[nodiscard]] auto headerData(int section, Qt::Orientation orientation, int role) const -> QVariant override;

private:
    QList<QString> headers_;
    QList<QList<QString>> rows_;
    QList<QString> tooltips_;
};

}  // namespace primerpicker::gui
