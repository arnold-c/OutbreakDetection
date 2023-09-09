#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using DifferentialEquations

include("ensemble-sim_single-scenario_infections.jl")

#%%
init_noise = [10.0]
sde_cb = DiscreteCallback(
    sde_condition, sde_affect!; save_positions = (false, false)
)

ensemble_noise_arr = create_static_noise_arr(
    init_noise,
    ensemble_single_scenario_spec.time_parameters,
    ensemble_single_scenario_spec.dynamics_parameters,
    ensemble_single_scenario_nsims;
    callback = sde_cb,
)
