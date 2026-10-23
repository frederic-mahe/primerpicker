// primerpicker-gui: entry point of the graphical front-end

#include "gui/application_icon.hpp"
#include "gui/main_window.hpp"

#include <QApplication>
#include <QCommandLineOption>
#include <QCommandLineParser>
#include <QMessageBox>
#include <QSettings>
#include <QString>

using namespace Qt::StringLiterals;

auto main(int argc, char* argv[]) -> int {
    auto const application = QApplication{argc, argv};
    QApplication::setOrganizationName(u"primerpicker"_s);
    QApplication::setApplicationName(u"primerpicker"_s);
    QApplication::setApplicationDisplayName(u"PrimerPicker"_s);
    QApplication::setApplicationVersion(QString::fromUtf8(PRIMERPICKER_VERSION));
    QApplication::setWindowIcon(primerpicker::gui::application_icon());

    auto parser = QCommandLineParser{};
    parser.setApplicationDescription(u"Interactive PCR and sequencing primer design"_s);
    parser.addHelpOption();
    parser.addVersionOption();
    auto const config = QCommandLineOption{
        u"config"_s,
        u"A settings file of Primer Designer 2.0 (pd2.cfg): its criteria, tolerances, primer length, search "
        u"ranges and weights replace the built-in defaults."_s,
        u"PD2.CFG"_s};
    parser.addOption(config);
    parser.process(application);

    auto settings = QSettings{};
    auto window = primerpicker::gui::MainWindow{settings};
    if (parser.isSet(config)) {
        if (auto const error = window.load_settings(parser.value(config)); error.has_value()) {
            QMessageBox::warning(&window, u"Settings"_s, *error + u"\nThe built-in defaults are used."_s);
        }
    }
    window.show();
    return QApplication::exec();
}
