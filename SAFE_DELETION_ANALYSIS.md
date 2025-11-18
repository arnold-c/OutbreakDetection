# Safe Deletion Analysis for OutbreakDetection Project

## Analysis Date
November 18, 2025

## Scope
Analysis of `scripts/optimal-thresholds.jl` and all its dependencies to identify unused functions and files.

## Scripts Analyzed
1. `scripts/optimal-thresholds.jl` (main entry point)
2. `scripts/plotting-setup.jl`
3. `scripts/schematic-plot.jl`
4. `scripts/optimal-thresholds_loading.jl`
5. `scripts/optimal-thresholds_plots.jl`
6. `scripts/optimal-thresholds_checks.jl`
7. `scripts/supplemental_tables.jl`
8. `scripts/supplemental_plots.jl`

## Files That Can Be Safely Removed

### 1. Source Files (src/)

#### src/plotting-functions.jl
- **Status**: SAFE TO DELETE
- **Reason**: File is essentially empty (only 1 line, likely just a newline)
- **Impact**: None - not used anywhere

#### src/makie-plotting-setup.jl
- **Status**: SAFE TO DELETE
- **Reason**: Not included or used in any of the analyzed scripts
- **Note**: Uses GLMakie instead of CairoMakie (scripts use CairoMakie)
- **Impact**: None for the manuscript workflow

### 2. Script Files (scripts/)

#### scripts/optimal-thresholds_tables.jl
- **Status**: SAFE TO DELETE
- **Reason**: File is empty (only 1 line)
- **Impact**: None - not actually doing anything

## Functions That Can Be Safely Removed from OutbreakDetectionUtils

### Functions NOT Used by the Scripts

Based on the analysis, the following functions are defined in `OutbreakDetectionUtils/src/` but are NOT used by any of the analyzed scripts:

#### From diag-testing-functions.jl
1. `calculate_daily_movingavg_startday`
2. `calculate_float_daily_movingavg`
3. `calculate_int_daily_movingavg`
4. `calculate_n_tests`
5. `calculate_period_sum`
6. `calculate_proportion_timeseries_in_outbreak` (the version taking outbreak_status_vec)
7. `calculate_outbreak_size`
8. `calculate_outbreak_duration` (both versions - the one taking lower/upper and the one taking row)
9. `match_outbreak_detection_bounds`
10. `calculate_delay_vec`
11. `freqtable_error_default_zero`

#### From optimal-threshold-functions.jl
1. `create_wide_optimal_thresholds_df`
2. `create_and_save_xlsx_optimal_threshold_summaries`
3. `save_xlsx_optimal_threshold_summaries`
4. `create_all_wide_optimal_threshold_summary_dfs`

#### From threshold-optimization-functions.jl
1. `setup_optimization`
2. `optimization_wrapper`

#### From scenario-optimizations.jl
1. `check_missing_scenario_optimizations`
2. `run_missing_scenario_optimizations!`
3. `get_most_recent_optimization_filepath`
4. `filter_optim_results`

#### From ensemble-functions.jl
1. `run_ensemble_jump_prob`
2. `run_jump_prob`
3. `summarize_ensemble_jump_prob`
4. `jump_prob_summary`
5. `run_OutbreakThresholdChars_creation`
6. `OutbreakThresholdChars_creation`
7. `get_ensemble_file` (all versions)
8. `match_ensemble_file!`
9. `collect_ensemble_file`

#### From cleaning-functions.jl
1. `create_sir_df_inner`
2. `create_sir_df` (both versions)
3. `create_sir_beta_dfs`
4. `create_sir_all_sim_quantiles`
5. `create_sir_all_sim_quantiles!`

#### From detection-thresholds.jl
1. `create_inc_infec_arr!`
2. `create_inc_infec_arr`

#### From noise-functions.jl
1. `create_noise_arr`

#### From collect-thresholds-vec_functions.jl
1. `collect_threshold_char_vec`

#### From SEIR-model.jl
1. `seir_mod_loop!`
2. `seir_mod!`
3. `convert_svec_to_array`
4. `convert_svec_to_matrix!`

#### From structs.jl
1. `ScenarioSpecification`
2. `TestPositivity`

#### From transmission-functions.jl
1. `calculate_beta` (the version taking R_0, gamma, mu, contact_mat, pop_matrix)
2. `calculateR0` (both versions)

## Functions USED by the Scripts

### Core Functions (DO NOT DELETE)

#### From OutbreakDetectionUtils:
- `StateParameters`
- `DynamicsParameters`
- `IndividualTestSpecification`
- `SimTimeParameters`
- `OutbreakSpecification`
- `OutbreakDetectionSpecification`
- `PoissonNoiseSpecification`
- `DynamicalNoiseSpecification`
- `MSO`
- `seir_mod`
- `convert_svec_to_matrix`
- `calculate_movingavg`
- `calculate_positives`
- `calculate_true_positives!`
- `calculate_noise_positives!`
- `detectoutbreak`
- `calculate_outbreak_thresholds`
- `calculate_outbreak_duration!`
- `classify_all_outbreaks!`
- `filter_only_outbreaks`
- `create_combinations_vec`
- `run_scenario_optimizations`
- `reshape_optim_df_to_matrix`
- `create_optimal_thresholds_df`
- `create_optimal_threshold_summary_df`
- `calculate_OptimalThresholdCharacteristics`
- `get_noise_description`
- `get_noise_magnitude`
- `get_test_description`
- `table_test_type`
- `plot_test_description`
- `arithmetic_mean`
- `calculate_f_beta_score`
- `calculate_tested!`
- `calculate_positives!`
- `calculate_noutbreaks`
- `calculate_outbreak_detection_characteristics`
- `calculate_OutbreakThresholdChars`
- `calculate_optimal_threshold`
- `calculate_ensemble_objective_metric`
- `calculate_outbreak_detection_accuracy`
- `objective_function`
- `run_optimization`

#### Constants (DO NOT DELETE):
- `BETA_MEAN`
- `BETA_FORCE`
- `SIGMA`
- `GAMMA`
- `MU`
- `ANNUAL_BIRTHS_PER_K`
- `EPSILON`
- `R0`
- `VACCINATION_COVERAGE`
- `POPULATION_N`
- `LATENT_PER_DAYS`
- `DUR_INF_DAYS`
- `LIFE_EXPECTANCY_YEARS`

#### From OutbreakDetection (src/):
- `line_plot`
- `collect_OptimalThresholdCharacteristics`
- Color constants (N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR, N_ALERTS_COLOR, etc.)

## Critical Issues Found

### Missing Variable Definitions
The following variables are used in `scripts/optimal-thresholds_loading.jl` but are NEVER DEFINED:
- `ensemble_state_specification`
- `ensemble_model_type`
- `ensemble_time_specification`
- `ensemble_nsims`
- `ensemble_noise_specification_vec`

**This is a bug** - these variables need to be defined before the script can run successfully.

## Recommendations

### Immediate Actions
1. **Delete empty/unused files**:
   - `src/plotting-functions.jl`
   - `src/makie-plotting-setup.jl`
   - `scripts/optimal-thresholds_tables.jl`

2. **Fix the missing variable definitions** in `scripts/optimal-thresholds_loading.jl` before attempting to run the script.

### Future Cleanup (Requires Testing)
The functions listed in "Functions NOT Used by the Scripts" can potentially be removed, but this requires:
1. Checking if they're used by other scripts not analyzed here
2. Checking if they're used by tests
3. Verifying they're not part of the public API that external users might depend on

### Conservative Approach
If you want to be very conservative:
1. Only delete the three empty/nearly-empty files listed above
2. Keep all functions in OutbreakDetectionUtils for now
3. Add deprecation warnings to unused functions before removing them

## Notes
- This analysis is based on static code analysis of the manuscript generation workflow
- Some functions may be used in tests or other workflows not analyzed here
- The `@static if false` block in `src/OutbreakDetection.jl` is for LSP support only and doesn't affect runtime
