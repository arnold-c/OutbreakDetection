#%%
using DrWatson
@quickactivate "OutbreakDetection"

include("ensemble-sim_single-scenario_noise.jl")

include(srcdir("makie-plotting-setup.jl"))

#%%
visualize_ensemble_noise(ensemble_noise_arr, time_p)
