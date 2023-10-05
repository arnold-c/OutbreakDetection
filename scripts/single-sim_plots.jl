#%%
using DrWatson
@quickactivate "OutbreakDetection"

using JLD2
using DataFrames
using DataFramesMeta

using OutbreakDetection

include(srcdir("makie-plotting-setup.jl"))

@unpack singlesim_time_p = load("data/singlesim/single-sim_setup.jld2")
@unpack trange = singlesim_time_p;

@unpack seir_array, inc_vec, beta_vec, seir_df = load(
    "data/singlesim/single-sim_arrays.jld2"
)

#%%
singlesim_timeseries_plot = draw_sir_plot(
    seir_df;
    annual = true,
    colors = seircolors,
    labels = ["S", "E", "I", "R", "N"]
)

save(plotsdir("singlesim/single-sim_timeseries.png"), singlesim_timeseries_plot)

#%%
singlesim_si_state_space_plot = @chain DataFrame(Tables.table(seir_array)) begin
    hcat(trange, _)
    rename!(["time", ["S", "E", "I", "R", "N"]...])
    data(_) *
    mapping(:I, :S; color = :time) *
    visual(Lines)
    draw
end

save(
    plotsdir("singlesim/single-sim_SI-state-space.png"),
    singlesim_si_state_space_plot,
)

#%%
singlesim_beta_plot = @chain DataFrame(Tables.table(beta_vec)) begin
    hcat(trange, _)
    rename!([:time, :beta_t])
    stack(_, Not("time"); variable_name = :beta, value_name = :Number)
    data(_) *
    mapping(
        :time => (t -> t / 365) => "Time (years)", :Number;
    ) *
    visual(Lines; linewidth = 1)
    draw(; facet = (; linkyaxes = :none), axis = (; limits = ((0, 3), nothing)))
end

save(plotsdir("singlesim/single-sim_beta.png"), singlesim_beta_plot)
