#pragma once

#include "core/sequence/nucleotide.hpp"
#include "core/sequence/topology.hpp"

#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QObject>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QString>
#include <QTabWidget>
#include <QVBoxLayout>
#include <QWidget>

#include <span>
#include <string_view>

namespace primerpicker::gui {

// Molecule-Ops -> View/Edit, read only for now: PRIMER.exe's View/Edit
// Molecule screens for the loaded molecule (tests/legacy/view_edit/):
// its header (name, size), its sequence in lines of 60 bases (Format
// is not recorded: single strand only), its features table and its
// annotations (the lines of its file before the sequence), and Find
// Sequence (IUPAC letters, both strands)
class MoleculeView final : public QWidget {
    Q_OBJECT

public:
    MoleculeView();

    struct Molecule {
        QString name;
        std::span<sequence::Nucleotide const> sequence;
        std::string_view file_text;  // the molecule file, for its features and annotations
        sequence::Topology topology;
    };
    void show_molecule(Molecule const& molecule);

    // Find Sequence: the sites of PATTERN (IUPAC letters), as listed
    void find(QString const& pattern);

private:
    std::span<sequence::Nucleotide const> molecule_;
    sequence::Topology topology_ = sequence::Topology::linear;
    // members are destroyed in reverse order: the widgets before the
    // layouts that hold them
    QVBoxLayout layout_;
    QFormLayout header_layout_;
    QHBoxLayout find_layout_;
    QLabel name_;
    QLabel size_;
    QTabWidget tabs_;
    QPlainTextEdit sequence_;
    QPlainTextEdit features_;
    QPlainTextEdit annotations_;
    QLabel search_label_;
    QLineEdit pattern_;
    QPushButton find_button_;
    QLabel sites_count_;
    QPlainTextEdit sites_;
};

}  // namespace primerpicker::gui
