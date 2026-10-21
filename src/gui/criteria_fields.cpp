#include "gui/criteria_fields.hpp"

#include "core/analysis/legacy_verdict.hpp"
#include "core/thermo/units.hpp"

#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QString>

#include <cstddef>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// the ranges of the fields
constexpr auto max_percent = 100;
constexpr auto min_celsius = -273;
constexpr auto max_celsius = 999;
constexpr auto max_threshold = 99;  // two digits
constexpr auto max_kcal = 99.9;
constexpr auto kcal_decimals = 1;
constexpr auto kcal_step = 0.1;

void set_up(QSpinBox& box, QString const& name, int const minimum, int const maximum) {
    box.setObjectName(name);
    box.setRange(minimum, maximum);
}

auto count(QSpinBox const& box) -> std::size_t {
    return static_cast<std::size_t>(box.value());
}

auto value(std::size_t const count) -> int {
    return static_cast<int>(count);
}

}  // namespace

CriteriaFields::CriteriaFields(QString const& prefix) {
    set_up(gc_min_, prefix + u"_gc_min_box"_s, 0, max_percent);
    set_up(gc_max_, prefix + u"_gc_max_box"_s, 0, max_percent);
    set_up(tm_min_, prefix + u"_tm_min_box"_s, min_celsius, max_celsius);
    set_up(tm_max_, prefix + u"_tm_max_box"_s, min_celsius, max_celsius);
    hairpin_.setObjectName(prefix + u"_hairpin_box"_s);
    hairpin_.setRange(-max_kcal, max_kcal);
    hairpin_.setDecimals(kcal_decimals);
    hairpin_.setSingleStep(kcal_step);
    hairpin_.setSuffix(u" kcal"_s);
    set_up(dimer_3prime_, prefix + u"_dimer_3prime_box"_s, 0, max_threshold);
    set_up(dimer_adjacent_, prefix + u"_dimer_adjacent_box"_s, 0, max_threshold);
    set_up(runs_, prefix + u"_runs_box"_s, 0, max_threshold);
    set_up(gc_run_, prefix + u"_gc_run_box"_s, 0, max_threshold);
    layout_.addRow(u"% GC: Min"_s, &gc_min_);
    layout_.addRow(u"% GC: Max"_s, &gc_max_);
    layout_.addRow(u"Tm °C: Min"_s, &tm_min_);
    layout_.addRow(u"Tm °C: Max"_s, &tm_max_);
    layout_.addRow(u"Hairpins: energy cutoff"_s, &hairpin_);
    layout_.addRow(u"Dimers: reject if >= matches at 3' end"_s, &dimer_3prime_);
    layout_.addRow(u"Dimers: reject if >= adjacent homologous bases"_s, &dimer_adjacent_);
    layout_.addRow(u"Runs: reject runs of >= bases (any)"_s, &runs_);
    layout_.addRow(u"Runs: reject runs of >= G or C at 3' end"_s, &gc_run_);
    setLayout(&layout_);
}

void CriteriaFields::set_criteria(analysis::LegacyCriteria const& criteria) {
    gc_min_.setValue(criteria.gc_min_percent);
    gc_max_.setValue(criteria.gc_max_percent);
    tm_min_.setValue(criteria.tm_min_celsius);
    tm_max_.setValue(criteria.tm_max_celsius);
    hairpin_.setValue(criteria.hairpin_cutoff.kcal_per_mol);
    dimer_3prime_.setValue(value(criteria.reject_dimer_3prime_matches));
    dimer_adjacent_.setValue(value(criteria.reject_dimer_adjacent_bases));
    runs_.setValue(value(criteria.reject_run_length));
    gc_run_.setValue(value(criteria.reject_gc_run_3prime));
}

auto CriteriaFields::criteria() const -> analysis::LegacyCriteria {
    return {.gc_min_percent = gc_min_.value(),
            .gc_max_percent = gc_max_.value(),
            .tm_min_celsius = tm_min_.value(),
            .tm_max_celsius = tm_max_.value(),
            .hairpin_cutoff = thermo::KilocaloriesPerMole{hairpin_.value()},
            .reject_dimer_3prime_matches = count(dimer_3prime_),
            .reject_dimer_adjacent_bases = count(dimer_adjacent_),
            .reject_run_length = count(runs_),
            .reject_gc_run_3prime = count(gc_run_)};
}

auto CriteriaFields::valid() const -> bool {
    return gc_min_.value() <= gc_max_.value() and tm_min_.value() <= tm_max_.value();
}

}  // namespace primerpicker::gui
