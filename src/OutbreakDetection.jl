module OutbreakDetection

using DrWatson
using GLMakie
using ColorSchemes
using UnPack
using DataFrames
using Chain
using NaNMath: NaNMath
using FLoops
using OutbreakDetectionUtils

include("plotting-helpers.jl")
export ACCURACY_COLOR, DAILY_SENSITIVITY_COLOR, DAILY_SPECIFICITY_COLOR,
    DAILY_PPV_COLOR, DAILY_NPV_COLOR, DETECTION_DELAY_COLOR,
    N_ALERTS_PER_OUTBREAK_COLOR, N_FALSE_ALERTS_COLOR, N_ALERTS_COLOR,
    N_OUTBREAKS_COLOR, N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR,
    PERC_OUTBREAKS_MISSED_COLOR, PERC_ALERTS_CORRECT_COLOR,
    PERC_ALERTS_FALSE_COLOR

include("single-sim_plots.jl")
export single_seir_plot, single_seir_statespace_plot, single_seir_beta_plot

include("quantile_plots.jl")
export create_seir_quantiles_plot

include("noise_plots.jl")
export visualize_ensemble_noise

include("ensemble-inspection_plots.jl")
export incidence_prevalence_plot,
    incidence_testing_plot, testing_plot, ensemble_outbreak_distribution_plot

include("outbreak-threshold-chars_plots.jl")
export ensemble_OTChars_plot, save_compare_ensemble_OTchars_plot,
    compare_ensemble_OTchars_plots, ensemble_outbreak_detect_diff_plot,
    test_positivity_distribution_plot

include("optimal-thresholds_plots.jl")
export compare_optimal_thresholds_chars_plot,
    create_optimal_thresholds_chars_plot,
    compare_optimal_thresholds_test_chars_plot,
    create_optimal_thresholds_test_chars_plot

include("single-scenario_plots.jl")
export singlescenario_test_positivity_plot,
    create_testing_related_plots, plot_all_single_scenarios

include("threshold_comparison_plots.jl")
export plot_all_threshold_comparisons

include("isocline_plots.jl")
export isocline_accuracy_plot

include("line_plots.jl")
export line_plot,
    collect_OptimalThresholdCharacteristics

@static if false
    # Include manuscript scripts for LSP support
    # These scripts are not loaded at runtime but help the LSP
    # recognize plotting functions and provide autocomplete
    include("../scripts/manuscript/optimal-thresholds.jl")
    include("../scripts/manuscript/optimal-thresholds_loading.jl")
    include("../scripts/manuscript/optimal-thresholds_plots.jl")
    include("../scripts/manuscript/optimal-thresholds_checks.jl")
    include("../scripts/manuscript/optimal-thresholds_tables.jl")
    include("../scripts/manuscript/supplemental_plots.jl")
    include("../scripts/manuscript/supplemental_tables.jl")
    include("../scripts/manuscript/plotting-setup.jl")
    include("../scripts/manuscript/schematic-plot.jl")
end

end
