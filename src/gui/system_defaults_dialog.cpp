#include "gui/system_defaults_dialog.hpp"

#include "core/io/pd2_cfg.hpp"
#include "core/sequence/legacy_direct_entry.hpp"

#include <QCheckBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QObject>
#include <QPushButton>
#include <QSpinBox>
#include <QString>
#include <QWidget>

#include <cstddef>

using namespace Qt::StringLiterals;

namespace primerpicker::gui {

namespace {

// PRIMER.exe's message for a refused value (from the strings of
// examples/PRIMER.exe)
constexpr auto illegal_value_message = u"Illegal value for defaults";

// the ranges of the fields: every value a pd2.cfg file may hold
constexpr auto min_primer_length = 2;
constexpr auto max_tolerance = 999;
constexpr auto max_offset = 9999;
constexpr auto max_weight = 5;

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

SystemDefaultsDialog::SystemDefaultsDialog(QWidget& parent)
    : QDialog{&parent},
      system_box_{u"System Defaults"_s},
      exclude_{u"False priming and uniqueness search: exclude Tm <= 0 °C"_s},
      sort_{u"Sort homologies by Tm"_s},
      criteria_box_{u"Criteria"_s},
      criteria_{u"system_defaults_dialog"_s},
      weights_box_{u"Find-Pair Weighting Factors (0 = ignore, 1-5)"_s},
      buttons_{QDialogButtonBox::Ok | QDialogButtonBox::Cancel},
      load_{u"&Load Settings File..."_s},
      save_{u"Save &As..."_s} {
    setObjectName(u"system_defaults_dialog"_s);
    setWindowTitle(u"System Defaults"_s);
    set_up(primer_length_, u"primer_length_box"_s, min_primer_length,
           static_cast<int>(sequence::direct_entry_max_length));
    exclude_.setObjectName(u"exclude_box"_s);
    sort_.setObjectName(u"sort_box"_s);
    set_up(gc_tolerance_, u"gc_tolerance_box"_s, 0, max_tolerance);
    set_up(tm_tolerance_, u"tm_tolerance_box"_s, 0, max_tolerance);
    set_up(match_start_, u"match_start_box"_s, 0, max_offset);
    set_up(match_end_, u"match_end_box"_s, 0, max_offset);
    set_up(find_pair_start_, u"find_pair_start_box"_s, 0, max_offset);
    set_up(find_pair_end_, u"find_pair_end_box"_s, 0, max_offset);
    set_up(weight_gc_, u"weight_gc_box"_s, 0, max_weight);
    set_up(weight_tm_, u"weight_tm_box"_s, 0, max_weight);
    set_up(weight_3prime_, u"weight_3prime_box"_s, 0, max_weight);
    set_up(weight_any_, u"weight_any_box"_s, 0, max_weight);
    set_up(weight_false_priming_, u"weight_false_priming_box"_s, 0, max_weight);
    message_.setObjectName(u"system_defaults_message"_s);
    buttons_.setObjectName(u"system_defaults_buttons"_s);
    load_.setObjectName(u"load_settings_button"_s);
    save_.setObjectName(u"save_settings_button"_s);

    system_layout_.addRow(u"Primer length (incl. 5' extension)"_s, &primer_length_);
    system_layout_.addRow(&exclude_);
    system_layout_.addRow(&sort_);
    system_layout_.addRow(u"Match primer GC % to within ±"_s, &gc_tolerance_);
    system_layout_.addRow(u"Match primer Tm to within ± (°C)"_s, &tm_tolerance_);
    system_layout_.addRow(u"Match Active Primer: start @ target -"_s, &match_start_);
    system_layout_.addRow(u"Match Active Primer: end @ target +"_s, &match_end_);
    system_layout_.addRow(u"Find Pair: (normal) start @ target -"_s, &find_pair_start_);
    system_layout_.addRow(u"Find Pair: (complement) end @ target +"_s, &find_pair_end_);
    system_box_.setLayout(&system_layout_);
    criteria_layout_.addWidget(&criteria_);
    criteria_box_.setLayout(&criteria_layout_);
    weights_layout_.addRow(u"GC mismatch"_s, &weight_gc_);
    weights_layout_.addRow(u"Tm mismatch"_s, &weight_tm_);
    weights_layout_.addRow(u"3' dimers"_s, &weight_3prime_);
    weights_layout_.addRow(u"Any dimers"_s, &weight_any_);
    weights_layout_.addRow(u"False priming"_s, &weight_false_priming_);
    weights_box_.setLayout(&weights_layout_);
    buttons_.addButton(&load_, QDialogButtonBox::ActionRole);
    buttons_.addButton(&save_, QDialogButtonBox::ActionRole);
    layout_.addWidget(&system_box_);
    layout_.addWidget(&criteria_box_);
    layout_.addWidget(&weights_box_);
    layout_.addWidget(&message_);
    layout_.addWidget(&buttons_);
    setLayout(&layout_);
    QObject::connect(&buttons_, &QDialogButtonBox::accepted, this, &SystemDefaultsDialog::check_and_accept);
    QObject::connect(&buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
    QObject::connect(&load_, &QPushButton::clicked, this, &SystemDefaultsDialog::load_requested);
    QObject::connect(&save_, &QPushButton::clicked, this, &SystemDefaultsDialog::save_requested);
}

void SystemDefaultsDialog::start(io::LegacySystemDefaults const& defaults) {
    started_ = defaults;
    primer_length_.setValue(value(defaults.primer_length));
    exclude_.setChecked(defaults.exclude_tm_at_most_zero);
    sort_.setChecked(defaults.sort_homologies_by_tm);
    gc_tolerance_.setValue(defaults.gc_tolerance_percent);
    tm_tolerance_.setValue(defaults.tm_tolerance_celsius);
    match_start_.setValue(value(defaults.match_start_offset));
    match_end_.setValue(value(defaults.match_end_offset));
    find_pair_start_.setValue(value(defaults.find_pair_start_offset));
    find_pair_end_.setValue(value(defaults.find_pair_end_offset));
    criteria_.set_criteria(defaults.criteria);
    weight_gc_.setValue(defaults.weights.gc_mismatch);
    weight_tm_.setValue(defaults.weights.tm_mismatch);
    weight_3prime_.setValue(defaults.weights.dimer_3prime);
    weight_any_.setValue(defaults.weights.dimer_any);
    weight_false_priming_.setValue(defaults.weights.false_priming);
    message_.clear();
    open();
}

auto SystemDefaultsDialog::defaults() const -> io::LegacySystemDefaults {
    auto defaults = started_;
    defaults.primer_length = count(primer_length_);
    defaults.exclude_tm_at_most_zero = exclude_.isChecked();
    defaults.sort_homologies_by_tm = sort_.isChecked();
    defaults.gc_tolerance_percent = gc_tolerance_.value();
    defaults.tm_tolerance_celsius = tm_tolerance_.value();
    defaults.match_start_offset = count(match_start_);
    defaults.match_end_offset = count(match_end_);
    defaults.find_pair_start_offset = count(find_pair_start_);
    defaults.find_pair_end_offset = count(find_pair_end_);
    defaults.criteria = criteria_.criteria();
    defaults.weights = {.gc_mismatch = weight_gc_.value(),
                        .tm_mismatch = weight_tm_.value(),
                        .dimer_3prime = weight_3prime_.value(),
                        .dimer_any = weight_any_.value(),
                        .false_priming = weight_false_priming_.value()};
    return defaults;
}

void SystemDefaultsDialog::check_and_accept() {
    if (not criteria_.valid()) {
        message_.setText(QString::fromUtf16(illegal_value_message));
        return;
    }
    accept();
    Q_EMIT defaults_entered();
}

}  // namespace primerpicker::gui
