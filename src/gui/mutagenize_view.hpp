#pragma once

#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"

#include <QDialog>
#include <QEvent>
#include <QHBoxLayout>
#include <QKeyEvent>
#include <QLabel>
#include <QListWidget>
#include <QObject>
#include <QPushButton>
#include <QString>
#include <QVBoxLayout>
#include <QWidget>

#include <cstddef>
#include <vector>

namespace primerpicker::gui {

// A list shown by a key of the Mutagenize screen (F1: the enzymes and
// their sites; F2: the codon table); typing a letter jumps to the first
// row that starts with it, Esc closes it
class MutagenizeListDialog final : public QDialog {
    Q_OBJECT

public:
    struct ListNames {
        QString title;
        QString object_name;
    };
    MutagenizeListDialog(QWidget& parent, ListNames const& names);

    // opens the dialog on ROWS, ROW being current (window-modal)
    void start(QList<QString> const& rows, int row);

private:
    QVBoxLayout layout_;
    QListWidget list_;
    QLabel help_;
};

// Modify -> Mutagenize: PRIMER.exe's Mutagenize Primer screen. The
// primer under its molecule bases, the translations of the molecule
// (O, R, I) and of the edited primer (M, U, T), and the values of the
// original and mutated primers. The arrows move the edit highlight, a
// typed base (A, C, G, T) replaces the one under it (the highlight
// stays); F1 lists the enzymes and their sites, F2 the codon table, I
// shows the analysis of the mutated primer, S saves it to the Primer
// List, Esc leaves the screen.
class MutagenizeView final : public QWidget {
    Q_OBJECT

public:
    MutagenizeView();

    struct Start {
        QString molecule_name;
        std::vector<sequence::Nucleotide> molecule;
        std::size_t position;  // the primer's 5' end (top strand)
        std::vector<sequence::Nucleotide> primer;
        analysis::LegacyCriteria criteria;
    };
    // Requires a primer that lies on the molecule at the position
    void start(Start start);

    [[nodiscard]] auto original() const -> std::vector<sequence::Nucleotide> const&;
    [[nodiscard]] auto mutated() const -> std::vector<sequence::Nucleotide> const&;
    [[nodiscard]] auto cursor() const -> std::size_t;  // 0-based

    void move_cursor(int bases);
    void type_base(sequence::Nucleotide base);
    // PRIMER.exe's "Primer Saved", after S
    void show_saved();

Q_SIGNALS:
    void save_requested();     // S
    void inspect_requested();  // I
    void exit_requested();     // Esc

protected:
    // F1 and F2 are the screen's keys here, not the window's shortcuts
    auto event(QEvent* event) -> bool override;
    void keyPressEvent(QKeyEvent* event) override;

private:
    void show_state();
    void show_enzymes();
    void show_codons();

    std::vector<sequence::Nucleotide> molecule_;
    std::size_t position_ = 1;
    std::vector<sequence::Nucleotide> original_;
    std::vector<sequence::Nucleotide> mutated_;
    std::size_t cursor_ = 0;
    analysis::LegacyCriteria criteria_;
    // members are destroyed in reverse order: the widgets before the
    // layouts that hold them
    QVBoxLayout layout_;
    QHBoxLayout buttons_layout_;
    QLabel header_;
    QLabel screen_;
    QLabel original_values_;
    QLabel mutated_values_;
    QLabel message_;
    QLabel keys_;
    QPushButton enzymes_button_;
    QPushButton codons_button_;
    QPushButton inspect_button_;
    QPushButton save_button_;
    MutagenizeListDialog enzymes_;
    MutagenizeListDialog codons_;
};

}  // namespace primerpicker::gui
