#%%
using DrWatson
@quickactivate "OutbreakDetection"

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(funsdir("plotting-functions.jl"))
includet(scriptsdir("ensemble-noise-sim.jl"))

#%%
visualize_ensemble_noise(ensemble_noise_arr, time_p)
