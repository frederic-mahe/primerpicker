#pragma once

#include <QDialog>
#include <QDialogButtonBox>
#include <QLabel>
#include <QObject>
#include <QPushButton>
#include <QTabWidget>
#include <QTextBrowser>
#include <QVBoxLayout>
#include <QWidget>

namespace primerpicker::gui {

// Help -> About PrimerPicker (not in PRIMER.exe): the version, the
// licence, the Qt version, and the licences of the packages
// (packaging/licenses/, built into the program): the third-party
// notices, the GNU GPL and the GNU LGPL, version 3
class AboutDialog final : public QDialog {
    Q_OBJECT

public:
    explicit AboutDialog(QWidget& parent);

private:
    // members are destroyed in reverse order: the widgets before the
    // tabs and the layout that hold them
    QVBoxLayout layout_;
    QTabWidget tabs_;
    QDialogButtonBox buttons_;
    QLabel about_;
    QTextBrowser third_party_;
    QTextBrowser gpl_;
    QTextBrowser lgpl_;
    QPushButton about_qt_;
};

}  // namespace primerpicker::gui
