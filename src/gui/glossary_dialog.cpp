#include "gui/glossary_dialog.hpp"

#include "core/glossary/glossary.hpp"
#include "core/glossary/glossary_entry.hpp"
#include "core/glossary/glossary_table.hpp"

#include <QDialog>
#include <QItemSelectionModel>
#include <QList>
#include <QModelIndex>
#include <QObject>
#include <QString>
#include <QWidget>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <string_view>
#include <utility>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

auto from_view(std::string_view const text) -> QString {
    return QString::fromUtf8(text.data(), static_cast<qsizetype>(text.size()));
}

auto titles() -> QList<QString> {
    auto list = QList<QString>{};
    for (auto const& entry : glossary::glossary_table) {
        list.append(from_view(entry.title));
    }
    return list;
}

}  // namespace

auto definition_text(std::string_view const name) -> QString {
    auto const entry = glossary::find_entry(name);
    return entry.has_value() ? from_view(glossary::definition(*entry)) : QString{};
}

GlossaryDialog::GlossaryDialog(QWidget& parent) : QDialog{&parent}, titles_{titles()} {
    setObjectName(u"glossary_dialog"_s);
    setWindowTitle(u"Glossary"_s);
    list_.setObjectName(u"glossary_list"_s);
    list_.setModel(&titles_);
    list_.setEditTriggers(QAbstractItemView::NoEditTriggers);
    text_.setObjectName(u"glossary_text"_s);
    layout_.addWidget(&list_, 1);
    layout_.addWidget(&text_, 3);
    setLayout(&layout_);
    constexpr auto width = 900;
    constexpr auto height = 600;
    resize(width, height);
    QObject::connect(list_.selectionModel(), &QItemSelectionModel::currentChanged, this,
                     [this](QModelIndex const& current) { show_row(current.row()); });
}

auto GlossaryDialog::show_entry(QString const& name) -> bool {
    auto const entry = glossary::find_entry(name.toStdString());
    auto row = 0;
    if (entry.has_value()) {
        auto const* const found = std::ranges::find_if(
            glossary::glossary_table, [&entry](glossary::Entry const& candidate) { return candidate.key == entry->key; });
        row = static_cast<int>(std::distance(glossary::glossary_table.begin(), found));
    }
    list_.setCurrentIndex(titles_.index(row));
    show_row(row);
    show();
    raise();
    return entry.has_value();
}

void GlossaryDialog::show_row(int const row) {
    if (row < 0 or std::cmp_greater_equal(row, glossary::glossary_table.size())) {
        return;
    }
    auto const& entry = glossary::glossary_table.at(static_cast<std::size_t>(row));
    text_.setMarkdown(u"# "_s + from_view(entry.title) + u"\n\n"_s + from_view(entry.text));
}

}  // namespace primerpicker::gui
