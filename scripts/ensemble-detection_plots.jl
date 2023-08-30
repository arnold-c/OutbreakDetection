#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using ColorSchemes
using GLMakie
using AlgebraOfGraphics

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(scriptsdir("ensemble-detection.jl"))
includet(funsdir("plotting-functions.jl"))

#%%
outbreakcols = [ColorSchemes.magma[i] for i in (200, 20)]

detect_outbreak_plot(
    inc_infec_arr,
    ensemble_seir_arr,
    param_dict[:time_p];
    colormap = outbreakcols,
    xlims = (90, 100),
    ylims_inc = (0, 150),
    ylims_periodsum = (0, 1000),
)
