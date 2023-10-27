"""
Note that each file is a module that defines the exported functions.
They are just listed here for convenience of sourcing one file.
"""
module OutbreakDetection

# using Reexport

include("transmission-functions.jl")
export calculate_beta,
    calculate_beta_amp, calculateR0, calculate_import_rate,
    calculate_mu
# @reexport using .TransmissionFunctions

include("structs.jl")
export SimTimeParameters, EnsembleSpecification, DynamicsParameters,
    StateParameters, OutbreakThresholdChars, OutbreakDetectionSpecification,
    OutbreakSpecification, IndividualTestSpecification, NoiseSpecification,
    ScenarioSpecification, TestPositivity
# @reexport using .ODStructs

include("SEIR-model.jl")
export seir_mod, seir_mod!, seir_mod_loop!,
    convert_svec_to_matrix, convert_svec_to_matrix!, convert_svec_to_array
# @reexport using .SEIRModel

include("cleaning-functions.jl")
export create_sir_df, create_sir_beta_dfs, create_sir_sim_array!,
    create_sir_all_sim_quantiles, create_sir_all_sim_quantiles!
# @reexport using .CleaningFunctions

include("bifurcation-functions.jl")
export birth_rate_bifurcation_simulation!, bifurcation_summary,
    beta_force_bifurcation_simulation!,
    birth_rate_beta_force_bifurcation_simulation!,
    birth_rate_beta_force_bifurcation_annual_summary,
    birth_rate_beta_force_bifurcation_cycle_summary

include("detection-thresholds.jl")
export create_inc_infec_arr,
    create_inc_infec_arr!,
    calculate_outbreak_thresholds
# @reexport using .DetectionThresholds

include("diag-testing-functions.jl")
export create_testing_arrs, create_testing_arrs!, calculate_tested!,
    calculate_positives!,
    calculate_true_positives!, calculate_noise_positives!,
    calculate_movingavg, calculate_movingavg!,
    detectoutbreak, detectoutbreak!, calculate_daily_detection_characteristics,
    calculate_noutbreaks, calculate_OutbreakThresholdChars,
    calculate_test_positivity, calculate_outbreak_detection_characteristics
# @reexport using .DiagTestingFunctions

include("ensemble-functions.jl")
export create_combinations_vec, create_ensemble_spec_combinations,
    run_ensemble_jump_prob, run_jump_prob,
    summarize_ensemble_jump_prob, jump_prob_summary,
    run_OutbreakThresholdChars_creation, OutbreakThresholdChars_creation,
    get_ensemble_file

# @reexport using .EnsembleFunctions

include("noise-functions.jl")
export create_poisson_noise_arr, create_poisson_noise_arr!
# @reexport using .NoiseFunctions

include("plotting-functions.jl")
export seircolors, seir_state_labels, create_sir_plot, draw_sir_plot,
    bifurcation_plot, bifurcation_heatmap,
    sir_quantiles_array_base_plot, create_sir_quantiles_plot, outbreakcols,
    detect_outbreak_plot, visualize_ensemble_noise, incidence_testing_plot,
    testing_plot, ensemble_outbreak_distribution_plot, ensemble_OTChars_plot,
    singlescenario_test_positivity_plot, test_positivity_distribution_plot,
    ensemble_outbreak_detect_diff_plot, compare_ensemble_OTchars_plots
# @reexport using .PlottingFunctions

end
