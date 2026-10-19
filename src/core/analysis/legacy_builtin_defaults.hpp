#pragma once

#include <cstddef>

// The built-in System Defaults of Primer Designer 2.0: the values PRIMER.exe
// shows on a fresh start, without a pd2.cfg file (tests/legacy/README.md,
// "Built-in defaults", recorded 2026-10-17: screens/system_defaults.png,
// screens/criteria_defaults.png, screens/weights_defaults.png). The default
// member values of the settings structures come from here.
namespace primerpicker::analysis::builtin {

// System Defaults
constexpr auto primer_length = std::size_t{20};  // bases, including the 5' extension
constexpr auto exclude_tm_at_most_zero = true;   // false priming and uniqueness search
constexpr auto sort_homologies_by_tm = true;
constexpr auto gc_tolerance_percent = 5;   // match primer GC % to within
constexpr auto tm_tolerance_celsius = 10;  // match primer Tm to within
constexpr auto match_start_offset = std::size_t{50};  // Match Active Primer: from target - offset
constexpr auto match_end_offset = std::size_t{50};    // to target + offset
constexpr auto find_pair_start_offset = std::size_t{100};  // Find Pair, primer A: from target - offset
constexpr auto find_pair_end_offset = std::size_t{100};    // primer B (complement): to target + offset

// Criteria
constexpr auto gc_min_percent = 50;
constexpr auto gc_max_percent = 60;
constexpr auto tm_min_celsius = 55;
constexpr auto tm_max_celsius = 80;
constexpr auto hairpin_cutoff_kcal_per_mol = 0.0;
constexpr auto reject_dimer_3prime_matches = std::size_t{3};  // reject if >=
constexpr auto reject_dimer_adjacent_bases = std::size_t{7};  // reject if >=
constexpr auto reject_run_length = std::size_t{3};            // reject if >=
constexpr auto reject_gc_run_3prime = std::size_t{3};         // reject if >=

// Find-Pair weighting factors (all five)
constexpr auto find_pair_weight = 1;

}  // namespace primerpicker::analysis::builtin
