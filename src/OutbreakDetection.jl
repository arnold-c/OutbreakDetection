"""
Note that each file is a module that defines the exported functions.
They are just listed here for convenience of sourcing one file.
"""
module OutbreakDetection

# using Reexport

include("transmission-functions.jl")
export calculate_beta, calculateR0, calculate_import_rate, calculate_mu
# @reexport using .TransmissionFunctions

include("structs.jl")
export SimTimeParameters, EnsembleSpecification, DynamicsParameters,
    StateParameters, OutbreakThresholdChars, OutbreakDetectionSpecification,
    OutbreakSpecification, IndividualTestSpecification, NoiseSpecification,
    ScenarioSpecification
# @reexport using .ODStructs

include("SEIR-model.jl")
export calculate_beta_amp, seir_mod, seir_mod!, seir_mod_loop!
# @reexport using .SEIRModel

include("cleaning-functions.jl")
export create_sir_df, create_sir_beta_dfs, create_sir_sim_array!,
    create_sir_all_sims_array, create_sir_all_sims_array!,
   create_sir_all_sim_quantiles, create_sir_all_sim_quantiles!
# @reexport using .CleaningFunctions

include("detection-thresholds.jl")
export create_inc_infec_arr, create_inc_infec_arr!, calculate_outbreak_thresholds
# @reexport using .DetectionThresholds

include("diag-testing-functions.jl")
export create_testing_arr, create_testing_arr!, calculate_tested!,
    calculate_pos, calculate_pos!, calculate_movingavg, calculate_movingavg!,
    detectoutbreak, detectoutbreak!, calculate_ot_characterstics,
    calculate_noutbreaks, calculate_OutbreakThresholdChars, create_combinations_vec,
    run_OutbreakThresholdChars_creation, OutbreakThresholdChars_creation,
    get_scenario_file
# @reexport using .DiagTestingFunctions

include("ensemble-functions.jl")
export run_ensemble_jump_prob, run_jump_prob, summarize_ensemble_jump_prob,
    jump_prob_summary, get_ensemble_file
# @reexport using .EnsembleFunctions

include("noise-functions.jl")
export create_noise_arr, create_noise_arr!, sde_affect!, sde_condition
# @reexport using .NoiseFunctions

include("plotting-functions.jl")
export seircolors, seir_state_labels, create_sir_plot, draw_sir_plot,
    sir_quantiles_array_base_plot, create_sir_quantiles_plot, outbreakcols,
    detect_outbreak_plot, visualize_ensemble_noise, incidence_testing_plot,
    testing_plot, ensemble_outbreak_distribution_plot, ensemble_OTChars_plot
# @reexport using .PlottingFunctions

end
