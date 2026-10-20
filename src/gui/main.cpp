// primerpicker-gui: entry point of the graphical front-end

#include "gui/main_window.hpp"

#include <QApplication>

auto main(int argc, char* argv[]) -> int {
    auto const application = QApplication{argc, argv};
    auto window = primerpicker::gui::MainWindow{};
    window.show();
    return QApplication::exec();
}
