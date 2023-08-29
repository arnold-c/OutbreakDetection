#%%
using DrWatson
@quickactivate "OutbreakDetection"

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(scriptsdir("ensemble-sim.jl"))
includet(funsdir("plotting-functions.jl"))
includet(scriptsdir("ensemble-sim.jl"))

#%%
create_sir_quantiles_plot(
    ensemble_seir_summary; labels = seir_state_labels, colors = seircolors,
    annual = true, caption = caption, tstep = param_dict[:tstep], xlims = (80, 100),
    ylims = (0, 1_000),
)

