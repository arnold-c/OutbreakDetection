module OutbreakDetection

using GLMakie
using ColorSchemes
using UnPack
using DataFrames
using Chain
using NaNMath: NaNMath
using FLoops

include("makie-plotting-setup.jl")

include("plotting-helpers.jl")

include("single-sim_plots.jl")
export single_seir_plot, single_seir_statespace_plot

include("plotting-functions.jl")
export seircolors,
    seir_state_labels, create_sir_plot, draw_sir_plot,
    ACCURACY_COLOR, DAILY_SENSITIVITY_COLOR, DAILY_SPECIFICITY_COLOR,
    DAILY_PPV_COLOR, DAILY_NPV_COLOR,
    DETECTION_DELAY_COLOR, N_ALERTS_PER_OUTBREAK_COLOR,
    N_FALSE_ALERTS_COLOR, N_ALERTS_COLOR, N_OUTBREAKS_COLOR,
    N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR,
    PERC_OUTBREAKS_MISSED_COLOR, PERC_ALERTS_CORRECT_COLOR,
    PERC_ALERTS_FALSE_COLOR,
    sir_quantiles_array_base_plot, create_sir_quantiles_plot,
    incidence_prevalence_plot, visualize_ensemble_noise, incidence_testing_plot,
    testing_plot, ensemble_outbreak_distribution_plot, ensemble_OTChars_plot,
    singlescenario_test_positivity_plot, test_positivity_distribution_plot,
    ensemble_outbreak_detect_diff_plot, save_compare_ensemble_OTchars_plot,
    compare_ensemble_OTchars_plots,
    compare_optimal_thresholds_chars_plot, create_optimal_thresholds_chars_plot,
    compare_optimal_thresholds_test_chars_plot,
    create_optimal_thresholds_test_chars_plot

include("ensemble-sim_single-scenario_plots.jl")
export create_testing_related_plots, plot_all_single_scenarios

include("threshold_comparison_plots.jl")
export collect_threshold_char_vec, plot_all_threshold_comparisons

end
