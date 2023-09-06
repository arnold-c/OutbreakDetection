#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using DifferentialEquations
using DataFrames

include("../src/OutbreakDetection.jl")
using .OutbreakDetection

include("ensemble-sim_single-scenario_infections.jl")

#%%
init_noise = [10.0]
sde_cb = DiscreteCallback(
    sde_condition, sde_affect!; save_positions = (false, false)
)

ensemble_noise_arr = create_noise_arr(
    ensemble_jump_arr,
    init_noise,
    time_p,
    ensemble_dynamics_p;
    callback = sde_cb,
)
