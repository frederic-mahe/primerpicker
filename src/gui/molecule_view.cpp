#include "gui/molecule_view.hpp"

#include "core/io/legacy_molecule_view.hpp"
#include "core/search/legacy_find_sequence.hpp"
#include "core/sequence/topology.hpp"

#include <QFontDatabase>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QObject>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QString>
#include <QWidget>

#include <string>
#include <vector>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

auto joined(std::vector<std::string> const& lines) -> QString {
    auto text = QList<QString>{};
    for (auto const& line : lines) {
        text.append(QString::fromStdString(line));
    }
    return text.join(u'\n');
}

void make_read_only(QPlainTextEdit& text, QString const& name) {
    text.setObjectName(name);
    text.setReadOnly(true);
    text.setLineWrapMode(QPlainTextEdit::NoWrap);
    text.setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
}

}  // namespace

MoleculeView::MoleculeView() : search_label_{u"Search for:"_s}, find_button_{u"&Find"_s} {
    setObjectName(u"molecule_view"_s);
    name_.setObjectName(u"molecule_name_label"_s);
    size_.setObjectName(u"molecule_size_label"_s);
    make_read_only(sequence_, u"molecule_sequence"_s);
    make_read_only(features_, u"molecule_features"_s);
    make_read_only(annotations_, u"molecule_annotations"_s);
    make_read_only(sites_, u"molecule_sites"_s);
    pattern_.setObjectName(u"find_sequence_field"_s);
    pattern_.setPlaceholderText(u"IUPAC letters: N any base; M A/C, R A/G, W A/T, S C/G, Y C/T, K G/T; "
                                u"B not A, D not C, H not G, V not T"_s);
    find_button_.setObjectName(u"find_sequence_button"_s);
    sites_count_.setObjectName(u"find_sequence_count"_s);
    header_layout_.addRow(u"Name:"_s, &name_);
    header_layout_.addRow(u"Size:"_s, &size_);
    tabs_.setObjectName(u"molecule_tabs"_s);
    tabs_.addTab(&sequence_, u"&Sequence"_s);
    tabs_.addTab(&features_, u"F&eatures"_s);
    tabs_.addTab(&annotations_, u"&Annotations"_s);
    find_layout_.addWidget(&search_label_);
    find_layout_.addWidget(&pattern_, 1);
    find_layout_.addWidget(&find_button_);
    find_layout_.addWidget(&sites_count_);
    QObject::connect(&find_button_, &QPushButton::clicked, this, [this] { find(pattern_.text()); });
    QObject::connect(&pattern_, &QLineEdit::returnPressed, this, [this] { find(pattern_.text()); });
    layout_.addLayout(&header_layout_);
    layout_.addWidget(&tabs_, 1);
    layout_.addLayout(&find_layout_);
    layout_.addWidget(&sites_);
    setLayout(&layout_);
}

void MoleculeView::show_molecule(Molecule const& molecule) {
    molecule_ = molecule.sequence;
    topology_ = molecule.topology;
    name_.setText(molecule.name);
    auto const topology = molecule.topology == sequence::Topology::circular ? u"circular"_s : u"linear"_s;
    size_.setText(u"%1 bp, DNA, %2"_s.arg(molecule.sequence.size()).arg(topology));
    sequence_.setPlainText(joined(io::legacy_sequence_lines(molecule.sequence)));
    features_.setPlainText(joined(io::legacy_feature_lines(molecule.file_text)));
    auto const annotations = io::legacy_annotation_lines(molecule.file_text);
    annotations_.setPlainText(annotations.empty() ? u"-- no annotations --"_s : joined(annotations));
    sites_count_.clear();
    sites_.clear();
    tabs_.setCurrentWidget(&sequence_);
}

void MoleculeView::find(QString const& pattern) {
    auto const sites = search::legacy_find_sequence(molecule_, pattern.toStdString(), topology_);
    if (not sites.has_value()) {
        sites_count_.setText(u"Not a sequence of IUPAC letters"_s);
        sites_.clear();
        return;
    }
    sites_count_.setText(u"Number of Sites = %1"_s.arg(sites->size()));  // PRIMER.exe's words
    sites_.setPlainText(joined(search::legacy_site_list_lines(*sites)));
}

}  // namespace primerpicker::gui
