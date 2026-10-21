#include "gui/main_window.hpp"
#include "gui_test_helpers.hpp"

#include <QAbstractButton>
#include <QAction>
#include <QApplication>
#include <QByteArray>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QIODevice>
#include <QMessageBox>
#include <QPushButton>
#include <QSpinBox>
#include <QString>
#include <QStringView>
#include <QTemporaryDir>
#include <QTimer>
#include <QWidget>

#include <catch2/catch_test_macros.hpp>

#include <functional>
#include <utility>

using namespace Qt::StringLiterals;

namespace {

using primerpicker::gui::MainWindow;
using primerpicker::gui::test::child;

auto legacy_file(QStringView const relative_path) -> QString {
    return QString::fromUtf8(PRIMERPICKER_LEGACY_DIR) + u'/' + relative_path;
}

auto read_all(QString const& path) -> QByteArray {
    auto file = QFile{path};
    REQUIRE(file.open(QIODevice::ReadOnly));
    return file.readAll();
}

// System Defaults: the primer length set to 25
void change_defaults(MainWindow& window) {
    child<QAction>(window, u"system_defaults").trigger();
    constexpr auto new_length = 25;
    child<QSpinBox>(window, u"primer_length_box").setValue(new_length);
    child<QDialogButtonBox>(window, u"system_defaults_buttons").button(QDialogButtonBox::Ok)->click();
}

// Runs ANSWER on the modal dialogs that open while the window closes
// (a message box, then a file dialog), until it returns true; it gives
// up after a while, closing any dialog, so that a test cannot hang
class DialogAnswers {
public:
    explicit DialogAnswers(std::function<bool(QWidget&)> answer) : answer_{std::move(answer)} {
        constexpr auto interval_ms = 20;
        timer_.setInterval(interval_ms);
        QObject::connect(&timer_, &QTimer::timeout, &timer_, [this] { poll(); });
        timer_.start();
    }

    [[nodiscard]] auto dialogs_seen() const -> int { return seen_; }

private:
    void poll() {
        constexpr auto give_up_after = 250;  // polls, 5 seconds
        auto* const dialog = QApplication::activeModalWidget();
        if (++polls_ > give_up_after) {
            if (dialog != nullptr) {
                dialog->close();
            }
            timer_.stop();
            return;
        }
        if (dialog == nullptr or dialog == last_) {
            return;
        }
        last_ = dialog;
        ++seen_;
        if (answer_(*dialog)) {
            timer_.stop();
        }
    }

    std::function<bool(QWidget&)> answer_;
    QTimer timer_;
    QWidget const* last_ = nullptr;
    int seen_ = 0;
    int polls_ = 0;
};

void click(QWidget& dialog, QMessageBox::StandardButton const button) {
    auto* const box = qobject_cast<QMessageBox*>(&dialog);
    REQUIRE(box != nullptr);
    CHECK(box->text() == u"System Defaults Changed But Not Saved!\nSave current defaults?"_s);
    box->button(button)->click();
}

}  // namespace

TEST_CASE("quitting without a change of the System Defaults asks nothing", "[gui][settings][quit]") {
    auto window = MainWindow{};
    window.show();
    CHECK_FALSE(window.settings_changed());
    auto const answers = DialogAnswers{[](QWidget& dialog) {
        dialog.close();
        return true;
    }};
    CHECK(window.close());
    CHECK(answers.dialogs_seen() == 0);
}

// Reference: PRIMER.exe asks "System Defaults Changed But Not Saved!
// SAVE CURRENT DEFAULTS (Y/N) ?" when it quits after a change
// (tests/legacy/README.md, "Settings file pd2.cfg")
TEST_CASE("after a change of the System Defaults, Cancel keeps the window open", "[gui][settings][quit]") {
    auto window = MainWindow{};
    window.show();
    change_defaults(window);
    CHECK(window.settings_changed());
    auto const answers = DialogAnswers{[](QWidget& dialog) {
        click(dialog, QMessageBox::Cancel);
        return true;
    }};
    CHECK_FALSE(window.close());
    CHECK(window.isVisible());
    CHECK(answers.dialogs_seen() == 1);
}

TEST_CASE("after a change of the System Defaults, No quits without saving", "[gui][settings][quit]") {
    auto window = MainWindow{};
    window.show();
    change_defaults(window);
    auto const answers = DialogAnswers{[](QWidget& dialog) {
        click(dialog, QMessageBox::No);
        return true;
    }};
    CHECK(window.close());
    CHECK(answers.dialogs_seen() == 1);
}

// Reference: the file PRIMER.exe writes when it quits after the primer
// length was set to 25 (tests/legacy/config/PD2_LENGTH25.CFG, recorded
// 2026-10-18)
TEST_CASE("after a change of the System Defaults, Yes saves them, then quits", "[gui][settings][quit][legacy]") {
    auto window = MainWindow{};
    window.show();
    change_defaults(window);
    auto const directory = QTemporaryDir{};
    REQUIRE(directory.isValid());
    auto const path = directory.filePath(u"PD2.CFG"_s);
    auto const answers = DialogAnswers{[&path](QWidget& dialog) {
        if (auto* const box = qobject_cast<QMessageBox*>(&dialog); box != nullptr) {
            click(dialog, QMessageBox::Yes);
            return false;  // the file dialog comes next
        }
        auto* const files = qobject_cast<QFileDialog*>(&dialog);
        REQUIRE(files != nullptr);
        files->selectFile(path);
        QDialog& file_dialog = *files;  // QFileDialog::accept is protected, QDialog's public
        file_dialog.accept();
        return true;
    }};
    CHECK(window.close());
    CHECK(answers.dialogs_seen() == 2);
    CHECK(read_all(path) == read_all(legacy_file(u"config/PD2_LENGTH25.CFG")));
}

TEST_CASE("changed System Defaults are saved as PD2.CFG, next to the settings file loaded", "[gui][settings][quit]") {
    auto window = MainWindow{};
    CHECK(window.settings_save_path() == QDir::current().filePath(u"PD2.CFG"_s));
    REQUIRE_FALSE(window.load_settings(legacy_file(u"config/PD2_CRITERIA.CFG")).has_value());
    CHECK_FALSE(window.settings_changed());
    CHECK(window.settings_save_path() == QDir{legacy_file(u"config"_s)}.filePath(u"PD2.CFG"_s));
}
