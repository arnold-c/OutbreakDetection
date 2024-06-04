module OutbreakDetection

using OutbreakDetectionUtils

export POPULATION_N, LATENT_PER_DAYS, DUR_INF_DAYS, R0, SIGMA, GAMMA,
    LIFE_EXPECTANCY_YEARS, ANNUAL_BIRTHS_PER_K, VACCINATION_COVERAGE,
    MU, BETA_MEAN, BETA_FORCE, EPSILON

export CLINICAL_CASE_TEST_SPEC, EPI_LINKED_CASE_TEST_SPEC, CLINICAL_TEST_SPECS

include("bifurcation-functions.jl")
export birth_rate_bifurcation_simulation!, bifurcation_summary,
    beta_force_bifurcation_simulation!,
    birth_rate_beta_force_bifurcation_simulation!,
    birth_rate_beta_force_bifurcation_annual_summary,
    birth_rate_beta_force_bifurcation_cycle_summary

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
    bifurcation_plot, bifurcation_heatmap,
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
