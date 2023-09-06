#%%
using DrWatson
@quickactivate "OutbreakDetection"

include("../src/OutbreakDetection.jl")
using .OutbreakDetection

include("ensemble-noise-sim.jl")

include(srcdir("makie-plotting-setup.jl"))

#%%
visualize_ensemble_noise(ensemble_noise_arr, time_p)
