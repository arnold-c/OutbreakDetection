#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack

include("ensemble-sim.jl")

#%%
ensemble_single_scenario_N = 500_000
ensemble_single_scenario_annual_births_per_k = 10
ensemble_single_scenario_beta_force = 0.2
ensemble_single_scenario_nsims = 1_000

#%%
ensemble_single_scenario_spec = EnsembleSpecification(
    ("seasonal-infectivity-import", "tau-leaping"),
    StateParameters(
        ensemble_single_scenario_N,
        Dict(
            :s_prop => 0.1,
            :e_prop => 0.01,
            :i_prop => 0.01,
            :r_prop => 0.88,
        ),
    ),
    DynamicsParameters(
        ensemble_single_scenario_N,
        ensemble_single_scenario_annual_births_per_k,
        ensemble_single_scenario_beta_force,
    ),
    SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0),
    ensemble_single_scenario_nsims
)

ensemble_single_scenario_sol = get_ensemble_file(
    "solution", ensemble_single_scenario_spec
)

ensemble_single_scenario_quantiles = get_ensemble_file(
    "95", ensemble_single_scenario_spec
)

@unpack ensemble_seir_arr, ensemble_jump_arr, ensemble_change_arr, dynamics_parameters, time_parameters, ensemble_param_dict = ensemble_single_scenario_sol

@unpack ensemble_seir_summary, caption = ensemble_single_scenario_quantiles

#%%
outbreak_spec = OutbreakSpecification(5, 30, 500)

inc_infec_arr = create_inc_infec_arr(ensemble_jump_arr, outbreak_spec)
