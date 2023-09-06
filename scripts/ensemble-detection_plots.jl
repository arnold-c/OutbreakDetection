#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ColorSchemes

include("../src/OutbreakDetection.jl")
using .OutbreakDetection

include("ensemble-detection.jl")

include(srcdir("makie-plotting-setup.jl"))

#%%
outbreakcols = [ColorSchemes.magma[i] for i in (200, 20)]

detect_outbreak_plot(
    inc_infec_arr,
    ensemble_seir_arr,
    ensemble_param_dict[:time_p];
    colormap = outbreakcols,
    xlims = (90, 100),
    ylims_inc = (0, 150),
    ylims_periodsum = (0, 1000),
)
