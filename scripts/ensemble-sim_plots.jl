#%%
using DrWatson
@quickactivate "OutbreakDetection"

include("../src/OutbreakDetection.jl")
using .OutbreakDetection

include("ensemble-sim.jl")

include(srcdir("makie-plotting-setup.jl"))

#%%
create_sir_quantiles_plot(
    ensemble_seir_summary; labels = seir_state_labels, colors = seircolors,
    annual = true, caption = caption, timeparams = time_p
)
