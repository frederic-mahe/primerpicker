#include "gui/about_dialog.hpp"

#include <QApplication>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFile>
#include <QIODevice>
#include <QObject>
#include <QString>
#include <QStringView>
#include <QTextBrowser>
#include <QTextEdit>
#include <QWidget>
#include <QtVersion>

#include <cassert>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// a file of the "licences" resources (src/gui/CMakeLists.txt)
auto resource_text(QString const& path) -> QString {
    auto file = QFile{path};
    [[maybe_unused]] auto const opened = file.open(QIODevice::ReadOnly);
    assert(opened);
    return QString::fromUtf8(file.readAll());
}

auto about_text() -> QString {
    return uR"(<h2>PrimerPicker %1</h2>
<p>Interactive PCR and sequencing primer design, by Frédéric Mahé.
Version 1 reproduces the features and results of Primer Designer 2.0
(Scientific &amp; Educational Software, 1990-1991); PrimerPicker is not
affiliated with Scientific &amp; Educational Software.</p>
<p>PrimerPicker is free software, under the GNU General Public License,
version 3 or later; it comes with no warranty.
<a href="https://github.com/frederic-mahe/primerpicker">github.com/frederic-mahe/primerpicker</a></p>
<p>Built with Qt %2, running with Qt %3 (GNU Lesser General Public
License, version 3: LGPL).</p>
<p>Restriction enzymes: REBASE (Roberts, Vincze, Posfai, Macelis, 2010.
Nucl. Acids Res. 38: D234-D236).</p>
<p>The licences of PrimerPicker and of the software it contains are in
the other tabs.</p>)"_s.arg(QApplication::applicationVersion(), QString::fromLatin1(QT_VERSION_STR),
                              QString::fromLatin1(qVersion()));
}

// a tab with the plain text of a licence
struct Licence {
    QStringView name;      // of the text browser
    QStringView resource;  // the file
};

void show_licence(QTextBrowser& browser, Licence const licence) {
    browser.setObjectName(licence.name.toString());
    browser.setLineWrapMode(QTextEdit::NoWrap);
    browser.setPlainText(resource_text(licence.resource.toString()));
}

}  // namespace

AboutDialog::AboutDialog(QWidget& parent) : QDialog{&parent} {
    setObjectName(u"about_dialog"_s);
    setWindowTitle(u"About PrimerPicker"_s);
    tabs_.setObjectName(u"about_tabs"_s);
    about_.setObjectName(u"about_text"_s);
    about_.setTextFormat(Qt::RichText);
    about_.setWordWrap(true);
    about_.setOpenExternalLinks(true);
    about_.setAlignment(Qt::AlignTop | Qt::AlignLeft);
    about_.setText(about_text());
    third_party_.setObjectName(u"third_party_text"_s);
    third_party_.setOpenExternalLinks(true);
    third_party_.setMarkdown(resource_text(u":/licences/packaging/licenses/THIRD-PARTY.md"_s));
    show_licence(gpl_, {.name = u"gpl_text", .resource = u":/licences/LICENSE"});
    show_licence(lgpl_, {.name = u"lgpl_text", .resource = u":/licences/packaging/licenses/LGPL-3.0.txt"});
    tabs_.addTab(&about_, u"PrimerPicker"_s);
    tabs_.addTab(&third_party_, u"Third-party software"_s);
    tabs_.addTab(&gpl_, u"GNU GPL 3"_s);
    tabs_.addTab(&lgpl_, u"GNU LGPL 3"_s);
    about_qt_.setObjectName(u"about_qt_button"_s);
    about_qt_.setText(u"About &Qt"_s);
    buttons_.addButton(&about_qt_, QDialogButtonBox::ActionRole);
    buttons_.addButton(QDialogButtonBox::Close);
    layout_.addWidget(&tabs_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    constexpr auto width = 720;
    constexpr auto height = 560;
    resize(width, height);
    QObject::connect(&about_qt_, &QPushButton::clicked, &QApplication::aboutQt);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
}

}  // namespace primerpicker::gui
