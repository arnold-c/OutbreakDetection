#%%
using DrWatson
@quickactivate "OutbreakDetection"

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(scriptsdir("ensemble-sim.jl"))
includet(funsdir("plotting-functions.jl"))
includet(scriptsdir("ensemble-sim.jl"))

GLMakie.activate!(; float = true)
set_aog_theme!()
update_theme!(resolution = (2200, 1300))

#%%
create_sir_quantiles_plot(
    ensemble_seir_summary; labels = seir_state_labels, colors = seircolors,
    annual = true, caption = caption, timeparams = time_p
)
