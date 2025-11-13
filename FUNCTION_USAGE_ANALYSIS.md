# Function Usage Analysis for src/ Directory

## Summary
This document provides a comprehensive analysis of all functions in the `src/` directory, identifying which are used by scripts and which are unused.

## Functions Used by Scripts

### line_plots.jl
- ✅ **line_plot** (2 methods) - USED
  - Called directly in: `scripts/manuscript/optimal-thresholds_plots.jl`, `scripts/manuscript/supplemental_plots.jl`
  - Usage: Main plotting function for manuscript figures

- ✅ **collect_OptimalThresholdCharacteristics** - USED
  - Called directly in: `scripts/manuscript/optimal-thresholds_loading.jl`
  - Called indirectly by: `line_plot` (first method)
  - Usage: Collects and structures optimal threshold data for plotting

- ✅ **_line_plot** - USED (indirectly)
  - Called by: `line_plot` (second method)
  - Usage: Internal helper for creating line plot facets

- ✅ **_line_plot_facet** - USED (indirectly)
  - Called by: `_line_plot`
  - Usage: Internal helper for individual facet creation

### plotting-helpers.jl
- ✅ **time_function** - POTENTIALLY USED
  - May be used by plotting functions

- ✅ **annual_label** - POTENTIALLY USED
  - May be used by plotting functions

- ✅ **lims_check** - POTENTIALLY USED
  - May be used by plotting functions

### Color Constants (plotting-helpers.jl)
- ✅ All color constants (ACCURACY_COLOR, etc.) - USED
  - Imported in: `scripts/manuscript/schematic-plot.jl`
  - Used throughout plotting functions

## Functions NOT Used by Scripts

### ensemble-inspection_plots.jl
- ❌ **incidence_prevalence_plot** - UNUSED
  - Exported but never called in scripts
  - Only called internally by `plot_all_single_scenarios`

- ❌ **incidence_testing_plot** - UNUSED
  - Exported but never called in scripts
  - Only called internally by `plot_all_single_scenarios`

- ❌ **testing_plot** - UNUSED
  - Exported but never called in scripts
  - Only called internally by `plot_all_single_scenarios`

- ❌ **ensemble_outbreak_distribution_plot** - UNUSED
  - Exported but never called in scripts
  - Only called internally by `plot_all_single_scenarios`

### noise_plots.jl
- ❌ **visualize_ensemble_noise** - UNUSED
  - Exported but never called in scripts
  - Only called internally by `plot_all_single_scenarios`

### single-sim_plots.jl
- ❌ **single_seir_plot** - UNUSED
  - Exported but never called in scripts or other src files

- ❌ **single_seir_statespace_plot** - UNUSED
  - Exported but never called in scripts or other src files

- ❌ **single_seir_beta_plot** - UNUSED
  - Exported but never called in scripts or other src files

### quantile_plots.jl
- ❌ **create_seir_quantiles_plot** - UNUSED
  - Exported but never called in scripts or other src files

- ❌ **seir_quantiles_array_base_plot** - UNUSED
  - Only called by `create_seir_quantiles_plot` (which is unused)

### outbreak-threshold-chars_plots.jl
- ❌ **ensemble_OTChars_plot** - UNUSED (by scripts)
  - Exported but never called directly in scripts
  - Only called internally by `plot_all_single_scenarios`

- ❌ **save_compare_ensemble_OTchars_plot** - UNUSED (by scripts)
  - Exported but never called directly in scripts
  - Only called internally by `plot_all_threshold_comparisons`

- ❌ **compare_ensemble_OTchars_plots** - UNUSED (by scripts)
  - Exported but never called directly in scripts
  - Only called internally by `save_compare_ensemble_OTchars_plot`

- ❌ **ensemble_outbreak_detect_diff_plot** - UNUSED (by scripts)
  - Exported but never called directly in scripts
  - Only called internally by `plot_all_single_scenarios`

- ❌ **test_positivity_distribution_plot** - UNUSED
  - Exported but never called in scripts or other src files

- ❌ **calculate_bins** - UNUSED (by scripts)
  - Only called internally by plotting functions

- ❌ **calculate_comparison_plot_facet_dims** - UNUSED (by scripts)
  - Only called internally by `compare_ensemble_OTchars_plots`

- ❌ **construct_single_OTchars_facet!** - UNUSED (by scripts)
  - Only called internally by plotting functions

- ❌ **construct_OTchars_facets!** - UNUSED (by scripts)
  - Only called internally by `compare_ensemble_OTchars_plots`

### optimal-thresholds_plots.jl
- ❌ **compare_optimal_thresholds_chars_plot** - UNUSED
  - Exported but never called in scripts or other src files

- ❌ **create_optimal_thresholds_chars_plot** - UNUSED (by scripts)
  - Exported but never called directly in scripts
  - Only called internally by `compare_optimal_thresholds_chars_plot`

- ❌ **compare_optimal_thresholds_test_chars_plot** - UNUSED
  - Exported but never called in scripts or other src files

- ❌ **create_optimal_thresholds_test_chars_plot** - UNUSED (by scripts)
  - Exported but never called directly in scripts
  - Only called internally by `compare_optimal_thresholds_test_chars_plot`

### single-scenario_plots.jl
- ❌ **singlescenario_test_positivity_plot** - UNUSED
  - Exported but never called in scripts or other src files

- ❌ **create_testing_related_plots** - UNUSED
  - Exported but never called in scripts or other src files

- ❌ **plot_all_single_scenarios** - UNUSED (by scripts)
  - Exported but never called directly in scripts
  - Only called internally by `create_testing_related_plots`

### threshold_comparison_plots.jl
- ❌ **plot_all_threshold_comparisons** - UNUSED
  - Exported but never called in scripts or other src files

### isocline_plots.jl
- ❌ **isocline_accuracy_plot** (2 methods) - UNUSED
  - Exported but never called in scripts or other src files

### ensemble-parameters.jl
- ❌ **All parameter definitions** - UNUSED (by scripts)
  - This file is included in `scripts/manuscript/optimal-thresholds_loading.jl`
  - But the parameters defined here may not be used

## Deletion Plan

Based on this analysis, the following files and functions can be safely deleted:

### Phase 1: Delete Completely Unused Files (No Script Dependencies)
These files contain functions that are never called by scripts, either directly or indirectly.

1. **single-sim_plots.jl** - All functions unused
2. **quantile_plots.jl** - All functions unused
3. **isocline_plots.jl** - All functions unused
4. **ensemble-parameters.jl** - Parameters not used by scripts

### Phase 2: Delete Unused High-Level Functions
These are exported functions that are never called by scripts.

5. **optimal-thresholds_plots.jl**:
   - Delete: `compare_optimal_thresholds_chars_plot`, `create_optimal_thresholds_chars_plot`
   - Delete: `compare_optimal_thresholds_test_chars_plot`, `create_optimal_thresholds_test_chars_plot`

6. **single-scenario_plots.jl**:
   - Delete: `singlescenario_test_positivity_plot`, `create_testing_related_plots`, `plot_all_single_scenarios`
   - Delete: All internal functions called only by these

7. **threshold_comparison_plots.jl**:
   - Delete: `plot_all_threshold_comparisons`

### Phase 3: Delete Unused Internal Functions
After deleting high-level functions, these internal functions become orphaned.

8. **ensemble-inspection_plots.jl**:
   - Delete: `incidence_prevalence_plot`, `incidence_testing_plot`, `testing_plot`, `ensemble_outbreak_distribution_plot`

9. **noise_plots.jl**:
   - Delete: `visualize_ensemble_noise`

10. **outbreak-threshold-chars_plots.jl**:
    - Delete: `ensemble_OTChars_plot`, `save_compare_ensemble_OTchars_plot`, `compare_ensemble_OTchars_plots`
    - Delete: `ensemble_outbreak_detect_diff_plot`, `test_positivity_distribution_plot`
    - Delete: `calculate_bins`, `calculate_comparison_plot_facet_dims`
    - Delete: `construct_single_OTchars_facet!`, `construct_OTchars_facets!`

### Files to Keep
- **line_plots.jl** - Core functionality used by manuscript scripts
- **plotting-helpers.jl** - Color constants and helper functions used throughout
- **OutbreakDetection.jl** - Module definition (will need updating after deletions)

## Recommended Deletion Order (JJ Revisions)

### Revision 1: Remove single-sim_plots.jl
**Type**: Structural change
**Files affected**:
- Delete: `src/single-sim_plots.jl`
- Edit: `src/OutbreakDetection.jl` (remove include and exports)

**Commit message**:
```
src: remove unused single simulation plotting functions

Delete single-sim_plots.jl containing single_seir_plot,
single_seir_statespace_plot, and single_seir_beta_plot. These functions
are exported but never called by any scripts or other source files.
```

### Revision 2: Remove quantile_plots.jl
**Type**: Structural change
**Files affected**:
- Delete: `src/quantile_plots.jl`
- Edit: `src/OutbreakDetection.jl` (remove include and exports)

**Commit message**:
```
src: remove unused quantile plotting functions

Delete quantile_plots.jl containing create_seir_quantiles_plot and
seir_quantiles_array_base_plot. These functions are exported but never
called by any scripts or other source files.
```

### Revision 3: Remove isocline_plots.jl
**Type**: Structural change
**Files affected**:
- Delete: `src/isocline_plots.jl`
- Edit: `src/OutbreakDetection.jl` (remove include and exports)

**Commit message**:
```
src: remove unused isocline plotting functions

Delete isocline_plots.jl containing two methods of isocline_accuracy_plot.
These functions are exported but never called by any scripts or other
source files.
```

### Revision 4: Remove ensemble-parameters.jl
**Type**: Structural change
**Files affected**:
- Delete: `src/ensemble-parameters.jl`
- Edit: `scripts/manuscript/optimal-thresholds_loading.jl` (remove include)

**Commit message**:
```
src: remove unused ensemble parameters file

Delete ensemble-parameters.jl as the parameters defined are not used by
manuscript scripts. The loading script includes this file but doesn't
reference the defined parameters.
```

### Revision 5: Remove unused optimal-thresholds_plots.jl functions
**Type**: Structural change
**Files affected**:
- Edit: `src/optimal-thresholds_plots.jl` (remove 4 functions)
- Edit: `src/OutbreakDetection.jl` (remove exports)

**Commit message**:
```
src: remove unused optimal threshold comparison functions

Remove compare_optimal_thresholds_chars_plot,
create_optimal_thresholds_chars_plot,
compare_optimal_thresholds_test_chars_plot, and
create_optimal_thresholds_test_chars_plot. These exported functions are
never called by scripts.
```

### Revision 6: Remove threshold_comparison_plots.jl
**Type**: Structural change
**Files affected**:
- Delete: `src/threshold_comparison_plots.jl`
- Edit: `src/OutbreakDetection.jl` (remove include and exports)

**Commit message**:
```
src: remove unused threshold comparison plotting module

Delete threshold_comparison_plots.jl containing
plot_all_threshold_comparisons. This function is exported but never
called by any scripts. It internally uses save_compare_ensemble_OTchars_plot
which will be removed in a subsequent revision.
```

### Revision 7: Remove single-scenario_plots.jl
**Type**: Structural change
**Files affected**:
- Delete: `src/single-scenario_plots.jl`
- Edit: `src/OutbreakDetection.jl` (remove include and exports)

**Commit message**:
```
src: remove unused single scenario plotting module

Delete single-scenario_plots.jl containing singlescenario_test_positivity_plot,
create_testing_related_plots, and plot_all_single_scenarios. These
exported functions are never called by scripts. This module internally
uses functions from ensemble-inspection_plots, noise_plots, and
outbreak-threshold-chars_plots which will be removed in subsequent
revisions.
```

### Revision 8: Remove ensemble-inspection_plots.jl
**Type**: Structural change
**Files affected**:
- Delete: `src/ensemble-inspection_plots.jl`
- Edit: `src/OutbreakDetection.jl` (remove include and exports)

**Commit message**:
```
src: remove unused ensemble inspection plotting module

Delete ensemble-inspection_plots.jl containing incidence_prevalence_plot,
incidence_testing_plot, testing_plot, and
ensemble_outbreak_distribution_plot. These functions were only called by
plot_all_single_scenarios which was removed in a previous revision.
```

### Revision 9: Remove noise_plots.jl
**Type**: Structural change
**Files affected**:
- Delete: `src/noise_plots.jl`
- Edit: `src/OutbreakDetection.jl` (remove include and exports)

**Commit message**:
```
src: remove unused noise visualization plotting module

Delete noise_plots.jl containing visualize_ensemble_noise. This function
was only called by plot_all_single_scenarios which was removed in a
previous revision.
```

### Revision 10: Remove outbreak-threshold-chars_plots.jl
**Type**: Structural change
**Files affected**:
- Delete: `src/outbreak-threshold-chars_plots.jl`
- Edit: `src/OutbreakDetection.jl` (remove include and exports)

**Commit message**:
```
src: remove unused outbreak threshold characteristics plotting module

Delete outbreak-threshold-chars_plots.jl containing ensemble_OTChars_plot,
save_compare_ensemble_OTchars_plot, compare_ensemble_OTchars_plots,
ensemble_outbreak_detect_diff_plot, test_positivity_distribution_plot,
and their helper functions. These functions were only called by
plot_all_single_scenarios and plot_all_threshold_comparisons which were
removed in previous revisions.
```

### Revision 11: Clean up optimal-thresholds_plots.jl (if empty)
**Type**: Structural change
**Files affected**:
- Delete: `src/optimal-thresholds_plots.jl` (if no functions remain)
- Edit: `src/OutbreakDetection.jl` (remove include if file deleted)

**Commit message**:
```
src: remove empty optimal-thresholds_plots.jl

Delete optimal-thresholds_plots.jl as all functions have been removed in
previous revisions.
```

## Final State

After all deletions, the `src/` directory will contain:

1. **OutbreakDetection.jl** - Module definition
2. **line_plots.jl** - Core line plotting functionality (used by manuscript)
3. **plotting-helpers.jl** - Color constants and helper functions

This represents a reduction from 14 files to 3 files, removing approximately 1,500+ lines of unused code.

## Testing Strategy

After each revision:
1. Run `just tests` to ensure no breakage in OutbreakDetectionUtils
2. Run manuscript scripts to ensure they still work:
   - `julia --project scripts/manuscript/optimal-thresholds.jl`
3. Check that the module still loads: `julia --project -e 'using OutbreakDetection'`

## Notes

- All deletions are structural changes (no behavior changes)
- Each revision should be tested independently before proceeding
- The order is designed to minimize dependencies (delete dependents before dependencies)
- Some functions like `calculate_bins` are internal helpers that will be removed when their parent functions are deleted
