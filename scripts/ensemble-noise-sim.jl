#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using DifferentialEquations
using DataFrames

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(funsdir("ensemble-functions.jl"))
includet(scriptsdir("ensemble-sim.jl"))
includet(funsdir("noise-functions.jl"))

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
