#pragma once

#include "core/analysis/legacy_verdict.hpp"

#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QObject>
#include <QSpinBox>
#include <QString>
#include <QWidget>

namespace primerpicker::gui {

// The fields of PRIMER.exe's Criteria form: GC % and Tm ranges, hairpin
// energy cutoff, the dimer and run thresholds ("reject if >="). Each
// field is named PREFIX_gc_min_box, PREFIX_gc_max_box, PREFIX_tm_min_box,
// PREFIX_tm_max_box, PREFIX_hairpin_box, PREFIX_dimer_3prime_box,
// PREFIX_dimer_adjacent_box, PREFIX_runs_box and PREFIX_gc_run_box. The
// ranges hold every value a pd2.cfg file may hold.
class CriteriaFields final : public QWidget {
    Q_OBJECT

public:
    explicit CriteriaFields(QString const& prefix);

    void set_criteria(analysis::LegacyCriteria const& criteria);
    [[nodiscard]] auto criteria() const -> analysis::LegacyCriteria;

    // whether each minimum is at most its maximum
    [[nodiscard]] auto valid() const -> bool;

private:
    // members are destroyed in reverse order: the fields before the
    // layout that holds them
    QFormLayout layout_;
    QSpinBox gc_min_;
    QSpinBox gc_max_;
    QSpinBox tm_min_;
    QSpinBox tm_max_;
    QDoubleSpinBox hairpin_;
    QSpinBox dimer_3prime_;
    QSpinBox dimer_adjacent_;
    QSpinBox runs_;
    QSpinBox gc_run_;
};

}  // namespace primerpicker::gui
