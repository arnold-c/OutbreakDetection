#%%
using DrWatson
@quickactivate "OutbreakDetection"
using Revise

include("ensemble-diag-testing_scenarios.jl")

includet(srcdir("makie-plotting-setup.jl"))

#%%

