#%%
using DrWatson
@quickactivate "OutbreakDetection"

include("../src/OutbreakDetection.jl")
using .OutbreakDetection

# Load infection time series
include("ensemble-sim_single-scenario_infections.jl")
include("ensemble-sim_single-scenario_noise.jl")


#%%
testlag = 3
moveavglag = 7
detectthreshold = 10
perc_clinic = 0.3
perc_clinic_test = 0.8
perc_tested = perc_clinic * perc_clinic_test
testsens = 0.9
testspec = 0.9

#%%
