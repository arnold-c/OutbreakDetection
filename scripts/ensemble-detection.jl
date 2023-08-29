#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(scriptsdir("ensemble-sim.jl"))
includet(funsdir("detection-thresholds.jl"))

#%%
outbreakthreshold = 5
minoutbreakdur = 30
minoutbreaksize = 500

inc_infec_arr = create_inc_infec_arr(
    ensemble_jump_arr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)
