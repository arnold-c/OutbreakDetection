function create_sir_plot(
    sol_df; colors=["dodgerblue4", "firebrick3", "chocolate2", "purple"]
)
    sir_plot =
        data(sol_df) *
        mapping(
            :time => "Time (days)", :Number; color=:State => sorter("S", "I", "R", "N")
        ) *
        visual(Lines; linewidth=4)

    return draw(sir_plot; palettes=(; color=colors))
end

function create_sir_quantiles_plot(sim_quantiles = sim_quantiles; lower, upper, quantiles, δt = δt)
    fig = Figure()
    ax = Axis(fig[1, 1])

    med_index = findfirst(isequal(0.5), quantiles)
    lower_index = findfirst(isequal(lower), quantiles)
    upper_index = findfirst(isequal(upper), quantiles)

    # Medians
    lines!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[med_index, :, 1];
        color=colors[1],
        linewidth=2,
        label="S",
    )
    lines!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[med_index, :, 2];
        color=colors[2],
        linewidth=2,
        label="I",
    )
    lines!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[med_index, :, 3];
        color=colors[3],
        linewidth=2,
        label="R",
    )
    lines!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[med_index, :, 4];
        color=colors[4],
        linewidth=2,
        label="N",
    )

    # User-specified quantiles
    band!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[lower_index, :, 1],
        sim_quantiles[upper_index, :, 1];
        color=(colors[1], 0.5),
    )
    band!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[lower_index, :, 2],
        sim_quantiles[upper_index, :, 2];
        color=(colors[2], 0.5),
    )
    band!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[lower_index, :, 3],
        sim_quantiles[upper_index, :, 3];
        color=(colors[3], 0.5),
    )
    band!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[lower_index, :, 4],
        sim_quantiles[upper_index, :, 4];
        color=(colors[4], 0.5),
    )

    Legend(fig[1, 2], ax, "State")

    return fig
end
