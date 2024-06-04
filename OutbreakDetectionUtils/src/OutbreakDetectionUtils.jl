module OutbreakDetectionUtils

include("DrWatson-helpers.jl")
export outdir


include("transmission-functions.jl")
export calculate_beta, calculate_beta_amp, calculateR0, calculate_import_rate, calculate_mu

include("structs.jl")
export SimTimeParameters, EnsembleSpecification, DynamicsParameters,
    StateParameters, OutbreakThresholdChars, OutbreakDetectionSpecification,
    OutbreakSpecification, IndividualTestSpecification,
    PoissonNoiseSpecification, DynamicalNoiseSpecification, NoiseSpecification,
    getdirpath,
    ScenarioSpecification, TestPositivity, OptimalThresholdCharacteristics

include("dynamics-constants.jl")
export POPULATION_N, LATENT_PER_DAYS, DUR_INF_DAYS, R0, SIGMA, GAMMA,
    LIFE_EXPECTANCY_YEARS, ANNUAL_BIRTHS_PER_K, VACCINATION_COVERAGE,
    MU, BETA_MEAN, BETA_FORCE, EPSILON

include("test-constants.jl")
export CLINICAL_CASE_TEST_SPEC, EPI_LINKED_CASE_TEST_SPEC, CLINICAL_TEST_SPECS


include("SEIR-model.jl")
export seir_mod, seir_mod!, seir_mod_loop!,
    convert_svec_to_matrix, convert_svec_to_matrix!, convert_svec_to_array

include("detection-thresholds.jl")
export create_inc_infec_arr, create_inc_infec_arr!, calculate_outbreak_thresholds

include("diag-testing-functions.jl")
export create_testing_arrs, create_testing_arrs!, calculate_tested!,
    calculate_positives!, calculate_true_positives!, calculate_noise_positives!,
    calculate_movingavg, calculate_movingavg!,
    detectoutbreak, detectoutbreak!, calculate_daily_detection_characteristics,
    calculate_noutbreaks, calculate_OutbreakThresholdChars,
    calculate_test_positivity, calculate_outbreak_detection_characteristics,
    filter_first_matched_bounds, calculate_first_matched_bounds_index,
    calculate_cases_before_after_alert!, calculate_cases_before_after_alert

include("ensemble-functions.jl")
export create_combinations_vec, create_ensemble_spec_combinations,
    run_ensemble_jump_prob, run_jump_prob,
    summarize_ensemble_jump_prob, jump_prob_summary,
    run_OutbreakThresholdChars_creation, OutbreakThresholdChars_creation,
    get_ensemble_file

include("noise-functions.jl")
export create_noise_arr, add_poisson_noise_arr!

include("optimal-threshold-functions.jl")
export calculate_optimal_threshold, calculate_OptimalThresholdCharacteristics,
    calculate_optimal_threshold_summaries,
    create_optimal_thresholds_df, create_wide_optimal_thresholds_df,
    create_and_save_xlsx_optimal_threshold_summaries,
    create_optimal_threshold_summary_df,
    create_wide_optimal_threshold_summary_df,
    create_all_wide_optimal_threshold_summary_dfs,
    save_xlsx_optimal_threshold_summaries,
    create_and_save_xlsx_optimal_threshold_summaries,
    gt_table

end
