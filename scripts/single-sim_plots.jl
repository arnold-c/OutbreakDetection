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

@unpack seir_array, change_array, jump_array, beta_arr, seir_df = load(
    "data/singlesim/single-sim_arrays.jld2"
)

#%%
singlesim_timeseries_plot = draw_sir_plot(
    seir_df;
    annual = true,
    colors = seircolors,
    labels = ["S", "E", "I", "R", "N", "incidence"]
)

save(plotsdir("singlesim/single-sim_timeseries.png"), singlesim_timeseries_plot)

#%%
singlesim_jump_plot = @chain DataFrame(Tables.table(jump_array)) begin
    hcat(trange, _)
    rename!([
        "time",
        "Infect",
        "Susceptible Births",
        "S death",
        "R death",
        "Import",
        "Protected Births",
        "Latent",
        "E death",
        "Recovery",
        "I death"
    ])
    stack(_, Not("time"); variable_name = :Jump, value_name = :Number)
    data(_) *
    mapping(
        :time => (t -> t / 365) => "Time (years)",
        :Number;
        color = :Jump,
        layout = :Jump,
    ) *
    visual(Lines; linewidth = 1)
    draw(;
        facet = (; linkyaxes = :none), axis = (; limits = (nothing, nothing))
    )
end

save(plotsdir("singlesim/single-sim_jumps.png"), singlesim_jump_plot)

#%%
change_labels = ["dS", "dE", "dI", "dR", "dN"]
singlesim_change_plot = @chain DataFrame(Tables.table(change_array)) begin
    hcat(trange, _)
    rename!([
        "time", change_labels...
    ])
    stack(_, Not("time"); variable_name = :Change, value_name = :Number)
    data(_) *
    mapping(
        :time => (t -> t / 365) => "Time (years)", :Number;
        color = :Change => sorter(change_labels...), layout = :Change,
    ) *
    visual(Lines; linewidth = 1)
    draw(;
        facet = (; linkyaxes = :none),
        palettes = (; color = seircolors),
        axis = (; limits = ((0, 100), nothing)),
    )
end

save(plotsdir("singlesim/single-sim_changes.png"), singlesim_change_plot)

#%%
singlesim_si_state_space_plot = @chain DataFrame(Tables.table(seir_array)) begin
    hcat(trange, _)
    rename!(["time", seir_state_labels...])
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
singlesim_beta_plot = @chain DataFrame(Tables.table(beta_arr)) begin
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
