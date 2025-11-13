# Detailed Deletion Plan for Unused src/ Functions

## Overview
This document provides a step-by-step plan for removing unused functions from the `src/` directory, with specific file changes and jj revision messages for each step.

## Pre-Deletion Verification

Before starting, verify current state:
```bash
# Ensure all tests pass
just tests

# Verify manuscript scripts work
julia --project scripts/manuscript/optimal-thresholds.jl

# Check module loads
julia --project -e 'using OutbreakDetection'
```

## Deletion Sequence

### Revision 1: Remove single-sim_plots.jl

**Rationale**: Contains 3 functions (single_seir_plot, single_seir_statespace_plot, single_seir_beta_plot) that are exported but never called anywhere.

**Files to modify**:
1. Delete: `src/single-sim_plots.jl`
2. Edit: `src/OutbreakDetection.jl`

**Changes to OutbreakDetection.jl**:
- Remove line 21-22:
  ```julia
  include("single-sim_plots.jl")
  export single_seir_plot, single_seir_statespace_plot, single_seir_beta_plot
  ```

**JJ revision message**:
```
src: remove unused single simulation plotting functions

Delete single-sim_plots.jl containing single_seir_plot,
single_seir_statespace_plot, and single_seir_beta_plot. These
functions are exported but never called by any scripts or other
source files.

This is a structural change with no behavioral impact.
```

**Verification**:
```bash
just tests
julia --project -e 'using OutbreakDetection'
```

---

### Revision 2: Remove quantile_plots.jl

**Rationale**: Contains 2 functions (create_seir_quantiles_plot, seir_quantiles_array_base_plot) that are exported but never called anywhere.

**Files to modify**:
1. Delete: `src/quantile_plots.jl`
2. Edit: `src/OutbreakDetection.jl`

**Changes to OutbreakDetection.jl**:
- Remove line 24-25:
  ```julia
  include("quantile_plots.jl")
  export create_seir_quantiles_plot
  ```

**JJ revision message**:
```
src: remove unused quantile plotting functions

Delete quantile_plots.jl containing create_seir_quantiles_plot and
seir_quantiles_array_base_plot. These functions are exported but
never called by any scripts or other source files.

This is a structural change with no behavioral impact.
```

**Verification**:
```bash
just tests
julia --project -e 'using OutbreakDetection'
```

---

### Revision 3: Remove isocline_plots.jl

**Rationale**: Contains 2 methods of isocline_accuracy_plot that are exported but never called anywhere.

**Files to modify**:
1. Delete: `src/isocline_plots.jl`
2. Edit: `src/OutbreakDetection.jl`

**Changes to OutbreakDetection.jl**:
- Remove line 52-53:
  ```julia
  include("isocline_plots.jl")
  export isocline_accuracy_plot
  ```

**JJ revision message**:
```
src: remove unused isocline plotting functions

Delete isocline_plots.jl containing two methods of
isocline_accuracy_plot. These functions are exported but never
called by any scripts or other source files.

This is a structural change with no behavioral impact.
```

**Verification**:
```bash
just tests
julia --project -e 'using OutbreakDetection'
```

---

### Revision 4: Remove ensemble-parameters.jl

**Rationale**: This file defines parameters but they are not used by the manuscript scripts. The file is included but its contents are not referenced.

**Files to modify**:
1. Delete: `src/ensemble-parameters.jl`
2. Edit: `scripts/manuscript/optimal-thresholds_loading.jl`

**Changes to optimal-thresholds_loading.jl**:
- Remove line 10:
  ```julia
  include(srcdir("ensemble-parameters.jl"))
  ```
- Remove lines 12-14 (the false block that references it):
  ```julia
  if false
      include("../src/ensemble-parameters.jl")
  end
  ```

**JJ revision message**:
```
src: remove unused ensemble parameters file

Delete ensemble-parameters.jl as the parameters defined are not
used by manuscript scripts. The loading script included this file
but didn't reference the defined parameters.

This is a structural change with no behavioral impact.
```

**Verification**:
```bash
just tests
julia --project scripts/manuscript/optimal-thresholds.jl
```

---

### Revision 5: Remove unused functions from optimal-thresholds_plots.jl

**Rationale**: Four exported functions are never called by scripts: compare_optimal_thresholds_chars_plot, create_optimal_thresholds_chars_plot, compare_optimal_thresholds_test_chars_plot, create_optimal_thresholds_test_chars_plot.

**Files to modify**:
1. Edit: `src/optimal-thresholds_plots.jl` (remove functions)
2. Edit: `src/OutbreakDetection.jl` (remove exports)

**Changes to optimal-thresholds_plots.jl**:
- Delete lines 3-63 (compare_optimal_thresholds_chars_plot function)
- Delete lines 65-149 (create_optimal_thresholds_chars_plot function)
- Delete lines 151-206 (compare_optimal_thresholds_test_chars_plot function)
- Delete lines 208-292 (create_optimal_thresholds_test_chars_plot function)
- Result: File should be empty or only contain imports

**Changes to OutbreakDetection.jl**:
- Remove lines 39-43:
  ```julia
  include("optimal-thresholds_plots.jl")
  export compare_optimal_thresholds_chars_plot,
      create_optimal_thresholds_chars_plot,
      compare_optimal_thresholds_test_chars_plot,
      create_optimal_thresholds_test_chars_plot
  ```

**JJ revision message**:
```
src: remove unused optimal threshold comparison functions

Remove compare_optimal_thresholds_chars_plot,
create_optimal_thresholds_chars_plot,
compare_optimal_thresholds_test_chars_plot, and
create_optimal_thresholds_test_chars_plot from
optimal-thresholds_plots.jl. These exported functions are never
called by scripts.

This is a structural change with no behavioral impact.
```

**Verification**:
```bash
just tests
julia --project scripts/manuscript/optimal-thresholds.jl
```

---

### Revision 6: Remove threshold_comparison_plots.jl

**Rationale**: Contains plot_all_threshold_comparisons which is exported but never called. It depends on save_compare_ensemble_OTchars_plot which will be removed later.

**Files to modify**:
1. Delete: `src/threshold_comparison_plots.jl`
2. Edit: `src/OutbreakDetection.jl`

**Changes to OutbreakDetection.jl**:
- Remove lines 49-50:
  ```julia
  include("threshold_comparison_plots.jl")
  export plot_all_threshold_comparisons
  ```

**JJ revision message**:
```
src: remove unused threshold comparison plotting module

Delete threshold_comparison_plots.jl containing
plot_all_threshold_comparisons. This function is exported but
never called by any scripts. It internally uses
save_compare_ensemble_OTchars_plot which will be removed in a
subsequent revision.

This is a structural change with no behavioral impact.
```

**Verification**:
```bash
just tests
julia --project -e 'using OutbreakDetection'
```

---

### Revision 7: Remove single-scenario_plots.jl

**Rationale**: Contains 3 exported functions (singlescenario_test_positivity_plot, create_testing_related_plots, plot_all_single_scenarios) that are never called by scripts. This module depends on functions from ensemble-inspection_plots, noise_plots, and outbreak-threshold-chars_plots.

**Files to modify**:
1. Delete: `src/single-scenario_plots.jl`
2. Edit: `src/OutbreakDetection.jl`

**Changes to OutbreakDetection.jl**:
- Remove lines 45-47:
  ```julia
  include("single-scenario_plots.jl")
  export singlescenario_test_positivity_plot,
      create_testing_related_plots, plot_all_single_scenarios
  ```

**JJ revision message**:
```
src: remove unused single scenario plotting module

Delete single-scenario_plots.jl containing
singlescenario_test_positivity_plot, create_testing_related_plots,
and plot_all_single_scenarios. These exported functions are never
called by scripts. This module internally uses functions from
ensemble-inspection_plots, noise_plots, and
outbreak-threshold-chars_plots which will be removed in
subsequent revisions.

This is a structural change with no behavioral impact.
```

**Verification**:
```bash
just tests
julia --project -e 'using OutbreakDetection'
```

---

### Revision 8: Remove ensemble-inspection_plots.jl

**Rationale**: Contains 4 functions (incidence_prevalence_plot, incidence_testing_plot, testing_plot, ensemble_outbreak_distribution_plot) that were only called by plot_all_single_scenarios, which was removed in Revision 7.

**Files to modify**:
1. Delete: `src/ensemble-inspection_plots.jl`
2. Edit: `src/OutbreakDetection.jl`

**Changes to OutbreakDetection.jl**:
- Remove lines 30-32:
  ```julia
  include("ensemble-inspection_plots.jl")
  export incidence_prevalence_plot,
      incidence_testing_plot, testing_plot, ensemble_outbreak_distribution_plot
  ```

**JJ revision message**:
```
src: remove unused ensemble inspection plotting module

Delete ensemble-inspection_plots.jl containing
incidence_prevalence_plot, incidence_testing_plot, testing_plot,
and ensemble_outbreak_distribution_plot. These functions were only
called by plot_all_single_scenarios which was removed in a
previous revision.

This is a structural change with no behavioral impact.
```

**Verification**:
```bash
just tests
julia --project -e 'using OutbreakDetection'
```

---

### Revision 9: Remove noise_plots.jl

**Rationale**: Contains visualize_ensemble_noise which was only called by plot_all_single_scenarios, which was removed in Revision 7.

**Files to modify**:
1. Delete: `src/noise_plots.jl`
2. Edit: `src/OutbreakDetection.jl`

**Changes to OutbreakDetection.jl**:
- Remove lines 27-28:
  ```julia
  include("noise_plots.jl")
  export visualize_ensemble_noise
  ```

**JJ revision message**:
```
src: remove unused noise visualization plotting module

Delete noise_plots.jl containing visualize_ensemble_noise. This
function was only called by plot_all_single_scenarios which was
removed in a previous revision.

This is a structural change with no behavioral impact.
```

**Verification**:
```bash
just tests
julia --project -e 'using OutbreakDetection'
```

---

### Revision 10: Remove outbreak-threshold-chars_plots.jl

**Rationale**: Contains 9 functions (ensemble_OTChars_plot, save_compare_ensemble_OTchars_plot, compare_ensemble_OTchars_plots, ensemble_outbreak_detect_diff_plot, test_positivity_distribution_plot, and helper functions) that were only called by plot_all_single_scenarios and plot_all_threshold_comparisons, both removed in previous revisions.

**Files to modify**:
1. Delete: `src/outbreak-threshold-chars_plots.jl`
2. Edit: `src/OutbreakDetection.jl`

**Changes to OutbreakDetection.jl**:
- Remove lines 34-37:
  ```julia
  include("outbreak-threshold-chars_plots.jl")
  export ensemble_OTChars_plot, save_compare_ensemble_OTchars_plot,
      compare_ensemble_OTchars_plots, ensemble_outbreak_detect_diff_plot,
      test_positivity_distribution_plot
  ```

**JJ revision message**:
```
src: remove unused outbreak threshold characteristics plotting

Delete outbreak-threshold-chars_plots.jl containing
ensemble_OTChars_plot, save_compare_ensemble_OTchars_plot,
compare_ensemble_OTchars_plots, ensemble_outbreak_detect_diff_plot,
test_positivity_distribution_plot, and their helper functions
(calculate_bins, calculate_comparison_plot_facet_dims,
construct_single_OTchars_facet!, construct_OTchars_facets!).
These functions were only called by plot_all_single_scenarios and
plot_all_threshold_comparisons which were removed in previous
revisions.

This is a structural change with no behavioral impact.
```

**Verification**:
```bash
just tests
julia --project -e 'using OutbreakDetection'
```

---

### Revision 11: Remove empty optimal-thresholds_plots.jl (if applicable)

**Rationale**: After Revision 5, this file may be empty or only contain imports. If so, it should be removed.

**Files to modify**:
1. Delete: `src/optimal-thresholds_plots.jl` (if empty)
2. Edit: `src/OutbreakDetection.jl` (if file deleted)

**Note**: Only proceed with this revision if the file is empty after Revision 5.

**JJ revision message**:
```
src: remove empty optimal-thresholds_plots.jl

Delete optimal-thresholds_plots.jl as all functions have been
removed in previous revisions, leaving only an empty file.

This is a structural change with no behavioral impact.
```

**Verification**:
```bash
just tests
julia --project -e 'using OutbreakDetection'
julia --project scripts/manuscript/optimal-thresholds.jl
```

---

## Final State

After all revisions, the `src/` directory structure will be:

```
src/
├── OutbreakDetection.jl          # Module definition (updated)
├── line_plots.jl                 # Core line plotting (kept)
└── plotting-helpers.jl           # Color constants and helpers (kept)
```

**Files deleted** (11 total):
1. single-sim_plots.jl
2. quantile_plots.jl
3. isocline_plots.jl
4. ensemble-parameters.jl
5. threshold_comparison_plots.jl
6. single-scenario_plots.jl
7. ensemble-inspection_plots.jl
8. noise_plots.jl
9. outbreak-threshold-chars_plots.jl
10. optimal-thresholds_plots.jl (if empty)
11. (potentially) optimal-thresholds_plots.jl

**Lines of code removed**: Approximately 1,500+ lines

## Post-Deletion Verification

After completing all revisions:

```bash
# Run full test suite
just tests

# Verify manuscript generation works
julia --project scripts/manuscript/optimal-thresholds.jl

# Check module loads correctly
julia --project -e 'using OutbreakDetection'

# Verify exports are correct
julia --project -e 'using OutbreakDetection; println(names(OutbreakDetection))'
```

Expected exports after all deletions:
- `line_plot`
- `collect_OptimalThresholdCharacteristics`
- Color constants (ACCURACY_COLOR, etc.)

## Rollback Strategy

If any revision causes issues:

```bash
# Undo the last revision
jj undo

# Or restore a specific revision
jj restore --from <revision-id>
```

## Notes

- Each revision is independent and can be tested separately
- All changes are structural (no behavior changes)
- The order minimizes dependencies (delete dependents before dependencies)
- Some helper functions (like `calculate_bins`) are removed with their parent functions
- The `plotting-helpers.jl` file is kept because it contains color constants used by the manuscript scripts
