// Runs the GUI tests: Catch2's session inside a QApplication (widgets
// need one). Without a display: the offscreen platform, unless
// QT_QPA_PLATFORM chooses another one (CTest also runs the program to
// list its tests, with the user's environment).

#include <QApplication>
#include <QByteArray>
#include <QtEnvironmentVariables>

#include <catch2/catch_session.hpp>

using namespace Qt::StringLiterals;

auto main(int argc, char* argv[]) -> int {
    constexpr auto platform_variable = "QT_QPA_PLATFORM";
    if (qEnvironmentVariableIsEmpty(platform_variable)) {
        qputenv(platform_variable, "offscreen"_ba);
    }
    auto const application = QApplication{argc, argv};
    return Catch::Session{}.run(argc, argv);
}
