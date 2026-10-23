#include "gui/return_key_filter.hpp"

#include <QEvent>
#include <QKeyEvent>
#include <QObject>

#include <cassert>
#include <functional>
#include <utility>

namespace primerpicker::gui {

ReturnKeyFilter::ReturnKeyFilter(std::function<void()> on_return) : on_return_{std::move(on_return)} {
    assert(on_return_);
}

auto ReturnKeyFilter::eventFilter(QObject* const watched, QEvent* const event) -> bool {
    if (event->type() != QEvent::KeyPress) {
        return QObject::eventFilter(watched, event);
    }
    auto const* const key_event = dynamic_cast<QKeyEvent const*>(event);
    if (key_event == nullptr || (key_event->key() != Qt::Key_Return && key_event->key() != Qt::Key_Enter)) {
        return QObject::eventFilter(watched, event);
    }
    on_return_();
    return true;
}

}  // namespace primerpicker::gui
