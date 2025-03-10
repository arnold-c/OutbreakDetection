#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using JLD2
using DataFrames
using DataFramesMeta

using OutbreakDetectionUtils
using OutbreakDetection
include(srcdir("makie-plotting-setup.jl"))

#%%
@unpack singlesim_time_p = load(
    joinpath(outdir("singlesim"), "single-sim_setup.jld2")
)
@unpack trange = singlesim_time_p;

@unpack seir_array, inc_vec, beta_vec, seir_df = load(
    joinpath(outdir("singlesim"), "single-sim_arrays.jld2")
)

#%%
singlesim_timeseries_plot = single_seir_plot(
    seir_df;
    annual = true,
)

mkpath(plotsdir("singlesim"))

save(plotsdir("singlesim/single-sim_timeseries.png"), singlesim_timeseries_plot)

#%%
singlesim_si_state_space_plot = single_seir_statespace_plot(seir_array)

save(
    plotsdir("singlesim/single-sim_SI-state-space.png"),
    singlesim_si_state_space_plot,
)

#%%
singlesim_beta_plot = single_seir_beta_plot(beta_vec, annual = false)

save(plotsdir("singlesim/single-sim_beta.png"), singlesim_beta_plot)
