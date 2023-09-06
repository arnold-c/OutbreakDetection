#%%
using DrWatson
@quickactivate "OutbreakDetection"

include("../src/OutbreakDetection.jl")
using .OutbreakDetection

include("ensemble-sim.jl")

#%%
outbreakthreshold = 5
minoutbreakdur = 30
minoutbreaksize = 500

inc_infec_arr = create_inc_infec_arr(
    ensemble_jump_arr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)
