#%%
using DrWatson
@quickactivate "OutbreakDetection"

includet(scriptsdir("single-sim.jl"))

#%%
draw_sir_plot(
    seir_df;
    annual = true,
    colors = seircolors,
    labels = seir_state_labels
)

#%%
@chain DataFrame(Tables.table(jump_array')) begin
    hcat(trange, _)
    rename!([
        "time", "Infect", "Latent", "Recov", "Birth", "S_death", "E_death",
        "I_death", "R_death",
        "Import",
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

#%%
change_labels = ["dS", "dE", "dI", "dR", "dN"]
@chain DataFrame(Tables.table(change_array')) begin
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

#%%
@chain DataFrame(Tables.table(seir_array')) begin
    hcat(trange, _)
    rename!(["time", seir_state_labels...])
    data(_) *
    mapping(:I, :S; color = :time) *
    visual(Lines)
    draw
end

#%%
@chain DataFrame(Tables.table(beta_arr)) begin
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

