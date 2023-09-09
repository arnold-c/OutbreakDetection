#%%
using DrWatson
@quickactivate "OutbreakDetection"

include("ensemble-sim_single-scenario.jl")

include(srcdir("makie-plotting-setup.jl"))

#%%
visualize_ensemble_noise(
    ensemble_single_scenario_spec.noise_specification.noise_array,
    ensemble_single_scenario_spec.noise_specification.time_parameters
)
