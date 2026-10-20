#pragma once

#include <QMainWindow>
#include <QObject>

namespace primerpicker::gui {

// The main window: Primer Designer's main menu as a menu bar. Items
// that are not implemented yet are disabled.
class MainWindow final : public QMainWindow {
    Q_OBJECT

public:
    MainWindow();
};

}  // namespace primerpicker::gui
