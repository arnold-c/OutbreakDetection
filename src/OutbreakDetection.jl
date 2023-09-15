"""
Note that each file is a module that defines the exported functions.
They are just listed here for convenience of sourcing one file.
"""
module OutbreakDetection

# using Reexport

include("transmission-functions.jl")
export calculate_beta, calculate_beta_amp, calculateR0, calculate_import_rate,
    calculate_mu
# @reexport using .TransmissionFunctions

include("structs.jl")
export SimTimeParameters, EnsembleSpecification, DynamicsParameters,
    StateParameters, OutbreakThresholdChars, OutbreakDetectionSpecification,
    OutbreakSpecification, IndividualTestSpecification, NoiseSpecification,
    ScenarioSpecification
# @reexport using .ODStructs

include("SEIR-model.jl")
export seir_mod, seir_mod!, seir_mod_loop!
# @reexport using .SEIRModel

include("cleaning-functions.jl")
export create_sir_df, create_sir_beta_dfs, create_sir_sim_array!,
    create_sir_all_sims_array, create_sir_all_sims_array!,
    create_sir_all_sim_quantiles, create_sir_all_sim_quantiles!
# @reexport using .CleaningFunctions

include("bifurcation-functions.jl")
export birth_rate_bifurcation_simulation!, bifurcation_summary,
    beta_force_bifurcation_simulation!

include("detection-thresholds.jl")
export create_inc_infec_arr,
    create_inc_infec_arr!, calculate_outbreak_thresholds,
    create_inc_infec_arr_long!
# @reexport using .DetectionThresholds

include("diag-testing-functions.jl")
export create_testing_arr, create_testing_arr!, calculate_tested!,
    calculate_pos, calculate_pos!, calculate_movingavg, calculate_movingavg!,
    detectoutbreak, detectoutbreak!, calculate_ot_characterstics,
    calculate_noutbreaks, calculate_OutbreakThresholdChars,
    run_OutbreakThresholdChars_creation, OutbreakThresholdChars_creation,
    get_scenario_file
# @reexport using .DiagTestingFunctions

include("ensemble-functions.jl")
export create_combinations_vec, create_ensemble_spec_combinations,
    run_ensemble_jump_prob, run_jump_prob, summarize_ensemble_jump_prob,
    jump_prob_summary, get_ensemble_file
# @reexport using .EnsembleFunctions

include("noise-functions.jl")
export create_static_noise_arr, create_static_noise_arr!, sde_affect!,
    sde_condition, create_static_NoiseSpecification
# @reexport using .NoiseFunctions

include("plotting-functions.jl")
export seircolors, seir_state_labels, create_sir_plot, draw_sir_plot,
    bifurcation_plot,
    sir_quantiles_array_base_plot, create_sir_quantiles_plot, outbreakcols,
    detect_outbreak_plot, visualize_ensemble_noise, incidence_testing_plot,
    testing_plot, ensemble_outbreak_distribution_plot, ensemble_OTChars_plot,
    compare_ensemble_OTchars_plots
# @reexport using .PlottingFunctions

end
