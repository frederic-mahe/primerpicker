#pragma once

#include <QEvent>
#include <QObject>

#include <functional>

namespace primerpicker::gui {

// Calls a function when Return or Enter is pressed in the widgets it
// filters, and stops the key there. Item views emit activated() for
// Return on Windows and Linux, but not on macOS (Return starts
// editing there): the filter gives every platform the same keys.
class ReturnKeyFilter final : public QObject {
public:
    explicit ReturnKeyFilter(std::function<void()> on_return);

protected:
    auto eventFilter(QObject* watched, QEvent* event) -> bool override;

private:
    std::function<void()> on_return_;
};

}  // namespace primerpicker::gui
