#pragma once

#include <QDialog>
#include <QHBoxLayout>
#include <QListView>
#include <QObject>
#include <QString>
#include <QStringListModel>
#include <QTextBrowser>
#include <QWidget>

#include <string_view>

namespace primerpicker::gui {

// The Definition paragraph of the built-in definition NAME (for
// tooltips); empty if there is no such entry
[[nodiscard]] auto definition_text(std::string_view name) -> QString;

// Help -> Glossary: the built-in definitions (data/glossary/glossary.md),
// a list of terms and the selected entry's text
class GlossaryDialog final : public QDialog {
    Q_OBJECT

public:
    explicit GlossaryDialog(QWidget& parent);

    // shows the dialog (not modal) at the entry NAME (a key or an alias),
    // or at the first entry; whether NAME was found
    auto show_entry(QString const& name) -> bool;

private:
    void show_row(int row);

    // members are destroyed in reverse order: the widgets before the
    // layout that holds them
    QHBoxLayout layout_;
    QStringListModel titles_;
    QListView list_;
    QTextBrowser text_;
};

}  // namespace primerpicker::gui
