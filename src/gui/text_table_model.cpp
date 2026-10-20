#include "gui/text_table_model.hpp"

#include <QAbstractTableModel>
#include <QList>
#include <QModelIndex>
#include <QString>
#include <QVariant>

#include <utility>

namespace primerpicker::gui {

TextTableModel::TextTableModel(QList<QString> headers) : headers_{std::move(headers)} {}

void TextTableModel::set_rows(QList<QList<QString>> rows) {
    beginResetModel();
    rows_ = std::move(rows);
    endResetModel();
}

auto TextTableModel::rowCount(QModelIndex const& parent) const -> int {
    return parent.isValid() ? 0 : static_cast<int>(rows_.size());
}

auto TextTableModel::columnCount(QModelIndex const& parent) const -> int {
    return parent.isValid() ? 0 : static_cast<int>(headers_.size());
}

auto TextTableModel::data(QModelIndex const& index, int const role) const -> QVariant {
    if (role != Qt::DisplayRole or not index.isValid() or index.row() >= rows_.size()) {
        return {};
    }
    auto const& row = rows_.at(index.row());
    return index.column() < row.size() ? QVariant{row.at(index.column())} : QVariant{};
}

auto TextTableModel::headerData(int const section, Qt::Orientation const orientation, int const role) const
    -> QVariant {
    if (role != Qt::DisplayRole or orientation != Qt::Horizontal or section < 0 or section >= headers_.size()) {
        return {};
    }
    return headers_.at(section);
}

}  // namespace primerpicker::gui
