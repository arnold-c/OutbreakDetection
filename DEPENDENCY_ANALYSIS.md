# OutbreakDetection Dependency Analysis

**Generated:** 2025-11-13  
**Cleanup Completed:** 2025-11-13

This document provides a comprehensive analysis of all functions and files used by the manuscript scripts, along with recommendations for files that can be safely deleted.

## Cleanup Status

✅ **Completed on 2025-11-13**

- Removed 16 unused scripts from `/scripts` directory
- Reorganized manuscript scripts from `manuscript/scripts/` to `scripts/manuscript/`
- Updated all path references in Justfile and manuscript scripts
- Updated LSP includes in `src/OutbreakDetection.jl` for manuscript scripts
- Created `manuscript/README.md` with generation instructions

All recommendations in this document have been implemented.

---

## Table of Contents

1. [Manuscript Script Flow](#manuscript-script-flow)
2. [Functions Used by Manuscript](#functions-used-by-manuscript)
3. [Functions Defined in Manuscript Scripts](#functions-defined-in-manuscript-scripts)
4. [Scripts Directory Analysis](#scripts-directory-analysis)
5. [Deletion Recommendations](#deletion-recommendations)

---

## Manuscript Script Flow

### Entry Point: `scripts/manuscript/optimal-thresholds.jl`

**Note:** Scripts were moved from `manuscript/scripts/` to `scripts/manuscript/` on 2025-11-13.

This is the main orchestrator that includes:

- [ ] **`plotting-setup.jl`** - Sets up CairoMakie theme and plotting parameters
- [ ] **`schematic-plot.jl`** - Creates the schematic figure
- [ ] **`optimal-thresholds_loading.jl`** - Loads and optimizes threshold data
- [ ] **`optimal-thresholds_plots.jl`** - Creates main manuscript plots
- [ ] **`supplemental_tables.jl`** - Creates supplemental tables
- [ ] **`supplemental_plots.jl`** - Creates supplemental plots
- [ ] **`optimal-thresholds_checks.jl`** - Validation checks

---

## Functions Used by Manuscript

### From `OutbreakDetectionUtils` Package

#### Core Optimization & Data Processing

- [ ] **`run_scenario_optimizations()`** - Main optimization runner
  - Used in: `optimal-thresholds_loading.jl:104`
  - Defined in: `OutbreakDetectionUtils/src/scenario-optimizations.jl`

- [ ] **`create_combinations_vec()`** - Creates parameter combinations
  - Used in: `optimal-thresholds_loading.jl:35, 48, 64, 78`
  - Defined in: `OutbreakDetectionUtils/src/ensemble-functions.jl`

- [ ] **`reshape_optim_df_to_matrix()`** - Reshapes optimization results
  - Used in: `optimal-thresholds_loading.jl:144, 146`
  - Defined in: `OutbreakDetectionUtils/src/scenario-optimizations.jl`

#### Struct Types

- [ ] **`IndividualTestSpecification`** - Test specifications
  - Used in: `optimal-thresholds_loading.jl:18, 23`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`DynamicsParameters`** - Dynamics parameters
  - Used in: `optimal-thresholds_loading.jl:36`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`EnsembleSpecification`** - Ensemble specifications
  - Used in: `optimal-thresholds_loading.jl:48`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`OutbreakSpecification`** - Outbreak definitions
  - Used in: `optimal-thresholds_loading.jl:64`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`OutbreakDetectionSpecification`** - Detection specifications
  - Used in: `optimal-thresholds_loading.jl:78`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`PoissonNoiseSpecification`** - Poisson noise specification
  - Used in: `optimal-thresholds_loading.jl:130`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`DynamicalNoiseSpecification`** - Dynamical noise specification
  - Used in: `optimal-thresholds_loading.jl:136`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`StateParameters`** - State parameters
  - Used in: `schematic-plot.jl:10, 27`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`SimTimeParameters`** - Time parameters
  - Used in: `schematic-plot.jl:46`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`OptimalThresholdCharacteristics`** - Optimal threshold characteristics
  - Used throughout plotting functions
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

#### Accuracy Functions

- [ ] **`arithmetic_mean()`** - Mean accuracy metric
  - Used in: `optimal-thresholds_loading.jl:101, 120`
  - Defined in: `OutbreakDetectionUtils/src/diag-testing-functions.jl`

- [ ] **`calculate_f_beta_score()`** - F-beta score metric
  - Used in: `optimal-thresholds_loading.jl:101, 125`
  - Defined in: `OutbreakDetectionUtils/src/diag-testing-functions.jl`

#### Data Frame Creation

- [ ] **`create_optimal_thresholds_df()`** - Creates threshold dataframe
  - Used in: `line_plots.jl:366` (called from supplemental_tables.jl)
  - Defined in: `OutbreakDetectionUtils/src/optimal-threshold-functions.jl`

- [ ] **`create_optimal_threshold_summary_df()`** - Creates summary dataframe
  - Used in: `supplemental_tables.jl:50, 59, 68`; `line_plots.jl:346`
  - Defined in: `OutbreakDetectionUtils/src/optimal-threshold-functions.jl`

#### Helper Functions

- [ ] **`get_test_description()`** - Gets test description string
  - Used in: `optimal-thresholds_checks.jl:18`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`get_noise_magnitude()`** - Gets noise magnitude
  - Used in: `optimal-thresholds_checks.jl:25`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`get_noise_description()`** - Gets noise description
  - Used in: `line_plots.jl:136, 156, 268`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`table_test_type()`** - Formats test type for tables
  - Used in: `supplemental_tables.jl:130`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

- [ ] **`plot_test_description()`** - Formats test for plots
  - Used in: `line_plots.jl:237`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

#### SEIR Model & Detection

- [ ] **`seir_mod()`** - SEIR model simulation
  - Used in: `schematic-plot.jl:127, 146`
  - Defined in: `OutbreakDetectionUtils/src/SEIR-model.jl`

- [ ] **`convert_svec_to_matrix()`** - Converts solution vector to matrix
  - Used in: `schematic-plot.jl:137, 154`
  - Defined in: `OutbreakDetectionUtils/src/SEIR-model.jl`

- [ ] **`calculate_movingavg()`** - Calculates moving average
  - Used in: `schematic-plot.jl:136, 153, 183`
  - Defined in: `OutbreakDetectionUtils/src/diag-testing-functions.jl`

- [ ] **`calculate_outbreak_thresholds()`** - Finds outbreak bounds
  - Used in: `schematic-plot.jl:72, 193`
  - Defined in: `OutbreakDetectionUtils/src/detection-thresholds.jl`

- [ ] **`classify_all_outbreaks!()`** - Classifies outbreaks
  - Used in: `schematic-plot.jl:76`
  - Defined in: `OutbreakDetectionUtils/src/detection-thresholds.jl`

- [ ] **`filter_only_outbreaks()`** - Filters to real outbreaks
  - Used in: `schematic-plot.jl:84`
  - Defined in: `OutbreakDetectionUtils/src/detection-thresholds.jl`

- [ ] **`calculate_positives()`** - Calculates test positives
  - Used in: `schematic-plot.jl:164, 172`
  - Defined in: `OutbreakDetectionUtils/src/detection-thresholds.jl`

- [ ] **`calculate_true_positives!()`** - True positives
  - Used in: `schematic-plot.jl:165`
  - Defined in: `OutbreakDetectionUtils/src/diag-testing-functions.jl`

- [ ] **`calculate_noise_positives!()`** - False positives
  - Used in: `schematic-plot.jl:173`
  - Defined in: `OutbreakDetectionUtils/src/diag-testing-functions.jl`

- [ ] **`detectoutbreak()`** - Detects outbreak alerts
  - Used in: `schematic-plot.jl:188`
  - Defined in: `OutbreakDetectionUtils/src/diag-testing-functions.jl`

- [ ] **`calculate_outbreak_duration!()`** - Calculates duration
  - Used in: `schematic-plot.jl:195`
  - Defined in: `OutbreakDetectionUtils/src/diag-testing-functions.jl`

#### Constants

- [ ] **`SIGMA`** - Latent period rate
  - Used in: `optimal-thresholds_loading.jl:41`
  - Defined in: `OutbreakDetectionUtils/src/dynamics-constants.jl`

- [ ] **`GAMMA`** - Recovery rate
  - Used in: `optimal-thresholds_loading.jl:42`
  - Defined in: `OutbreakDetectionUtils/src/dynamics-constants.jl`

- [ ] **`MSO`** - Optimization method
  - Used in: `optimal-thresholds_loading.jl:110`
  - Defined in: `OutbreakDetectionUtils/src/structs.jl`

### From `OutbreakDetection` Package

#### Plotting Functions

- [ ] **`line_plot()`** - Main line plot function
  - Used in: `optimal-thresholds_plots.jl:2, 27, 52, 76`
  - Used in: `supplemental_plots.jl:2, 26, 51, 76, 100, 123, 146, 171, 195`
  - Defined in: `src/line_plots.jl`

- [ ] **`collect_OptimalThresholdCharacteristics()`** - Collects characteristics for plotting
  - Used in: `line_plots.jl:49, 261`
  - Defined in: `src/line_plots.jl`

#### Color Constants

- [ ] **`N_MISSED_OUTBREAKS_COLOR`** - Color for missed outbreaks
  - Used in: `schematic-plot.jl:4-5`
  - Defined in: `src/plotting-helpers.jl`

- [ ] **`PERC_OUTBREAKS_DETECTED_COLOR`** - Color for detected outbreaks
  - Used in: `schematic-plot.jl:4-5`
  - Defined in: `src/plotting-helpers.jl`

- [ ] **`N_ALERTS_COLOR`** - Color for alerts
  - Used in: `schematic-plot.jl:4-5`
  - Defined in: `src/plotting-helpers.jl`

---

## Functions Defined in Manuscript Scripts

### `supplemental_tables.jl`

- [ ] **`create_wide_df()`** (line 91)
  - Purpose: Converts long to wide format for tables
  - Parameters: `long_df`, `outcome::Symbol`, `noise_order`, `digits`
  - Returns: Wide format DataFrame

### `optimal-thresholds_checks.jl`

- [ ] **`compare_optimal_solution_mean_extrema()`** (lines 5, 49, 72, 99)
  - Purpose: Multiple dispatch versions for validation
  - Four different method signatures for different input types

- [ ] **`extract_test_optimal_solutions()`** (line 119)
  - Purpose: Extracts test-specific solutions
  - Parameters: `optimal_solutions`, `test_spec`
  - Returns: Filtered outbreak threshold characteristics

### `schematic-plot.jl`

- [ ] **`get_outbreak_status()`** (line 63)
  - Purpose: Determines outbreak status from incidence vector
  - Parameters: `inc_vec`, `outbreak_specification`
  - Returns: `outbreak_status`, `outbreak_bounds`

- [ ] **`shift_vec()`** (line 98)
  - Purpose: Shifts vector by offset
  - Parameters: `invec`, `shift::Integer`
  - Returns: Shifted vector

- [ ] **`create_schematic_simulation()`** (line 117)
  - Purpose: Creates simulation for schematic figure
  - Parameters: Multiple state, dynamics, and specification parameters
  - Returns: Tuple of simulation results

- [ ] **`plot_schematic()`** (line 222)
  - Purpose: Plots schematic figure
  - Parameters: Simulation results and plotting options
  - Returns: Makie Figure

---

## Scripts Directory Analysis

### Files Referenced in Code

#### Actively Used (via `@static if false` in OutbreakDetection.jl)

These are currently disabled but referenced:

- [ ] **`scripts/debugging.jl`**
  - Referenced in: `src/OutbreakDetection.jl:63`
  - Status: Disabled via `@static if false`

- [ ] **`scripts/ensemble-diag-testing_scenarios_plots.jl`**
  - Referenced in: `src/OutbreakDetection.jl:61`
  - Status: Disabled via `@static if false`

- [ ] **`scripts/ensemble-sim_optimal-accuracy-lineplot.jl`**
  - Referenced in: `src/OutbreakDetection.jl:62`
  - Status: Disabled via `@static if false`

- [ ] **`scripts/single-sim_plots.jl`**
  - Referenced in: `src/OutbreakDetection.jl:60`
  - Status: Disabled via `@static if false`

#### Referenced by Other Scripts

- [ ] **`scripts/single-sim_setup.jl`**
  - Referenced in: `scripts/single-sim.jl`, `scripts/single-sim_plots.jl`
  - Status: Used by other scripts (which are themselves not used)

- [ ] **`scripts/threshold-optimization.jl`**
  - Referenced in: `OutbreakDetectionUtils/src/scenario-optimizations.jl` (in comments)
  - Status: May be documentation/example

### Files NOT Referenced Anywhere

These scripts appear to be standalone exploratory/development scripts:

- [ ] **`scripts/calculate-dynamical-noise-vaccintion-rates.jl`**
  - Purpose: Likely calculates vaccination rates for dynamical noise
  - References: None found

- [ ] **`scripts/ensemble-diag-testing_constant-thresholds.jl`**
  - Purpose: Ensemble testing with constant thresholds
  - References: None found

- [ ] **`scripts/ensemble-diag-testing_optimal-thresholds_single-timeseries.jl`**
  - Purpose: Optimal thresholds for single time series
  - References: None found

- [ ] **`scripts/ensemble-diag-testing_optimal-thresholds.jl`**
  - Purpose: Optimal thresholds for ensemble
  - References: None found

- [ ] **`scripts/ensemble-sim_noise-visualizations.jl`**
  - Purpose: Visualizing noise in ensemble simulations
  - References: None found

- [ ] **`scripts/ensemble-sim_optimal-accuracy-isocline.jl`**
  - Purpose: Isocline plots for optimal accuracy
  - References: None found

- [ ] **`scripts/ensemble-sim_single-scenario.jl`**
  - Purpose: Single scenario ensemble simulation
  - References: None found

- [ ] **`scripts/ensemble-sim.jl`**
  - Purpose: General ensemble simulation
  - References: None found

- [ ] **`scripts/line-plots.jl`**
  - Purpose: Line plotting (superseded by `src/line_plots.jl`?)
  - References: None found

- [ ] **`scripts/outbreak-detection-schematic.jl`**
  - Purpose: Outbreak detection schematic (superseded by manuscript version?)
  - References: None found

- [ ] **`scripts/schematic-plots.jl`**
  - Purpose: Schematic plots (superseded by manuscript version?)
  - References: None found

- [ ] **`scripts/single-sim.jl`**
  - Purpose: Single simulation runner
  - References: None found

---

## Deletion Recommendations

### High Confidence - Safe to Delete

These files are NOT referenced anywhere in the codebase and appear to be exploratory scripts:

- [ ] `scripts/calculate-dynamical-noise-vaccintion-rates.jl`
- [ ] `scripts/ensemble-diag-testing_constant-thresholds.jl`
- [ ] `scripts/ensemble-diag-testing_optimal-thresholds_single-timeseries.jl`
- [ ] `scripts/ensemble-diag-testing_optimal-thresholds.jl`
- [ ] `scripts/ensemble-sim_noise-visualizations.jl`
- [ ] `scripts/ensemble-sim_optimal-accuracy-isocline.jl`
- [ ] `scripts/ensemble-sim_single-scenario.jl`
- [ ] `scripts/ensemble-sim.jl`
- [ ] `scripts/line-plots.jl` (functionality moved to `src/line_plots.jl`)
- [ ] `scripts/outbreak-detection-schematic.jl` (superseded by `manuscript/scripts/schematic-plot.jl`)
- [ ] `scripts/schematic-plots.jl` (superseded by `manuscript/scripts/schematic-plot.jl`)
- [ ] `scripts/single-sim.jl`

### Medium Confidence - Potentially Deletable

These are currently disabled via `@static if false` but may have been useful for development:

- [ ] `scripts/debugging.jl` - May be useful for future debugging
- [ ] `scripts/ensemble-diag-testing_scenarios_plots.jl` - May be useful for exploratory analysis
- [ ] `scripts/ensemble-sim_optimal-accuracy-lineplot.jl` - May be useful for exploratory analysis
- [ ] `scripts/single-sim_plots.jl` - May be useful for exploratory analysis
- [ ] `scripts/single-sim_setup.jl` - Used by single-sim_plots.jl

### Keep - Referenced or Potentially Useful

- [ ] `scripts/threshold-optimization.jl` - Referenced in comments, may serve as documentation/example

---

## Critical Dependencies - MUST KEEP

### All OutbreakDetectionUtils Source Files

All files in `OutbreakDetectionUtils/src/` are actively used:

- [ ] `scenario-optimizations.jl` - Core optimization logic
- [ ] `optimal-threshold-functions.jl` - Threshold calculations
- [ ] `structs.jl` - All struct definitions
- [ ] `diag-testing-functions.jl` - Testing and detection
- [ ] `SEIR-model.jl` - Disease model
- [ ] `detection-thresholds.jl` - Outbreak detection
- [ ] `ensemble-functions.jl` - Ensemble utilities
- [ ] `noise-functions.jl` - Noise generation
- [ ] `threshold-optimization-functions.jl` - Optimization functions
- [ ] `transmission-functions.jl` - Transmission calculations
- [ ] `cleaning-functions.jl` - Data cleaning
- [ ] `collect-thresholds-vec_functions.jl` - Threshold collection
- [ ] `dynamics-constants.jl` - Constants
- [ ] `test-constants.jl` - Test constants
- [ ] `DrWatson-helpers.jl` - DrWatson utilities
- [ ] `OutbreakDetectionUtils.jl` - Main module file

### All OutbreakDetection Source Files

All files in `src/` are actively used:

- [ ] `line_plots.jl` - Used by manuscript
- [ ] `plotting-helpers.jl` - Color constants
- [ ] `ensemble-parameters.jl` - Ensemble parameters
- [ ] `OutbreakDetection.jl` - Main module file
- [ ] All other plotting files (even if not directly used by manuscript, they're part of the package API)

### All Manuscript Scripts

All files in `manuscript/scripts/` are actively used:

- [ ] `optimal-thresholds.jl` - Main entry point
- [ ] `plotting-setup.jl` - Theme setup
- [ ] `schematic-plot.jl` - Schematic figure
- [ ] `optimal-thresholds_loading.jl` - Data loading
- [ ] `optimal-thresholds_plots.jl` - Main plots
- [ ] `supplemental_tables.jl` - Tables
- [ ] `supplemental_plots.jl` - Supplemental plots
- [ ] `optimal-thresholds_checks.jl` - Validation

---

## Next Steps

1. [ ] Review each file marked for deletion individually
2. [ ] Check git history to understand when/why each file was created
3. [ ] Consider archiving rather than deleting (move to `scripts/archive/`)
4. [ ] Run manuscript build after any deletions to ensure nothing breaks
5. [ ] Update documentation to reflect any changes

---

## Notes

- This analysis was performed on 2025-11-13
- All line numbers are approximate and may shift with code changes
- Before deleting any files, consider creating a git branch for safety
- Some scripts may have historical value even if not currently used
