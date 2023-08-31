#%%
using DrWatson
@quickactivate "OutbreakDetection"

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(scriptsdir("ensemble-diag-testing-sim.jl"))
includet(funsdir("plotting-functions.jl"))

#%%
incidence_testing_plot(
    inc_infec_arr,
    testing_arr,
    time_p,
    detectthreshold;
    sim = 1,
)
