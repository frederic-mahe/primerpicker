#include "gui/multi_primer_view.hpp"

#include "core/analysis/legacy_multi_primer.hpp"
#include "core/analysis/legacy_verdict.hpp"
#include "core/sequence/nucleotide.hpp"
#include "gui/text_table_model.hpp"

#include <QAbstractItemView>
#include <QHeaderView>
#include <QList>
#include <QString>
#include <QTableView>

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

constexpr auto field_names = std::array{u'A', u'B', u'C', u'D'};

auto letters(std::span<sequence::Nucleotide const> const nucleotides) -> QString {
    auto text = QString{};
    for (auto const nucleotide : nucleotides) {
        text.append(QChar::fromLatin1(sequence::to_letter(nucleotide)));
    }
    return text;
}

template <typename Value>
auto flagged(Value const value, bool const flag) -> QString {
    return QString::number(value) + (flag ? u"!"_s : QString{});
}

auto matrix_headers() -> QList<QString> {
    auto headers = QList<QString>{QString{}};
    for (auto const name : field_names) {
        headers.append(QString{name});
    }
    return headers;
}

void set_up(QTableView& table, TextTableModel& model, QString const& name) {
    table.setObjectName(name);
    table.setModel(&model);
    table.setEditTriggers(QAbstractItemView::NoEditTriggers);
    table.verticalHeader()->hide();
    table.horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
}

enum struct DimerValue : std::uint8_t { three_prime, any };

// a cell of a matrix: empty for an empty field, a dot below the
// diagonal, else the value with its flag
auto cell_text(analysis::LegacyMultiPrimer const& analysis, std::size_t const row, std::size_t const column,
               DimerValue const value) -> QString {
    auto const both = analysis.columns.at(row).has_value() and analysis.columns.at(column).has_value();
    auto const& cell = analysis.dimers.at(row).at(column);
    if (both and column < row) {
        return u"."_s;
    }
    if (not both or not cell.has_value()) {
        return {};
    }
    return value == DimerValue::three_prime ? flagged(cell->dimer_3prime_matches, cell->flags.three_prime)
                                            : flagged(cell->dimer_adjacent_bases, cell->flags.adjacent);
}

// a matrix: a row per field (empty fields too), the field's name first
auto matrix_rows(analysis::LegacyMultiPrimer const& analysis, DimerValue const value) -> QList<QList<QString>> {
    auto rows = QList<QList<QString>>{};
    for (auto row = std::size_t{0}; row < analysis::legacy_multi_primer_fields; ++row) {
        auto texts = QList<QString>{QString{field_names.at(row)}};
        for (auto column = std::size_t{0}; column < analysis::legacy_multi_primer_fields; ++column) {
            texts.append(cell_text(analysis, row, column, value));
        }
        rows.append(texts);
    }
    return rows;
}

}  // namespace

MultiPrimerView::MultiPrimerView()
    : balance_title_{u"Primer Balance and Cautions"_s},
      balance_model_{{u"Field"_s, u"Primer"_s, u"Length"_s, u"% GC"_s, u"Tm °C"_s, u"Runs of G/C at 3' end"_s,
                      u"Runs of bases (any)"_s, u"Hairpins"_s}},
      three_prime_title_{u"Primer Dimers: 3' end"_s},
      three_prime_model_{matrix_headers()},
      any_title_{u"Primer Dimers: Any"_s},
      any_model_{matrix_headers()} {
    setObjectName(u"multi_primer_view"_s);
    set_up(balance_table_, balance_model_, u"multi_balance_table"_s);
    set_up(three_prime_table_, three_prime_model_, u"multi_3prime_table"_s);
    set_up(any_table_, any_model_, u"multi_any_table"_s);
    layout_.addWidget(&balance_title_);
    layout_.addWidget(&balance_table_);
    layout_.addWidget(&three_prime_title_);
    layout_.addWidget(&three_prime_table_);
    layout_.addWidget(&any_title_);
    layout_.addWidget(&any_table_);
    setLayout(&layout_);
}

void MultiPrimerView::show_primers(analysis::MultiPrimerFields const& fields, analysis::LegacyCriteria const& criteria) {
    auto const analysis = analysis::legacy_multi_primer(fields, criteria);
    auto balance = QList<QList<QString>>{};
    for (auto index = std::size_t{0}; index < analysis::legacy_multi_primer_fields; ++index) {
        auto const& column = analysis.columns.at(index);
        auto const& field = fields.at(index);
        if (not column.has_value() or not field.has_value()) {
            continue;
        }
        auto const& actual = column->actual;
        auto const& cautions = column->cautions;
        balance.append({QString{field_names.at(index)}, letters(*field), QString::number(column->length),
                        flagged(actual.gc_percent, cautions.gc), flagged(actual.tm_celsius, cautions.tm),
                        flagged(actual.gc_run_3prime, cautions.gc_run_3prime),
                        flagged(actual.longest_run, cautions.runs), cautions.hairpin ? u"!"_s : u"OK"_s});
    }
    balance_model_.set_rows(std::move(balance));
    three_prime_model_.set_rows(matrix_rows(analysis, DimerValue::three_prime));
    any_model_.set_rows(matrix_rows(analysis, DimerValue::any));
}

}  // namespace primerpicker::gui
