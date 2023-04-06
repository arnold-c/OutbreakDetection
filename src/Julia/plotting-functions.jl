function create_sir_plot(sol_df)
    return data(sol_df) *
           mapping(
               :time => "Time (days)", :Number; color = :State => sorter("S", "I", "R", "N")
           ) *
           visual(Lines; linewidth = 4)
end

function create_annual_sir_plot(sol_df)
    return data(sol_df) *
           mapping(
               :time => (t -> t / 365.0) => "Time (years)", :Number;
               color = :State => sorter("S", "I", "R", "N"),
           ) *
           visual(Lines; linewidth = 4)
end

function draw_sir_plot(
    sir_plot; colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]
)
    return draw(sir_plot; palettes = (; color = colors))
end

function draw_sir_plot(
    sol_df::DataFrame;
    colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"],
    annual = false,
)
    return if annual == false
        draw_sir_plot(create_sir_plot(sol_df); colors = colors)
    else
        draw_sir_plot(create_annual_sir_plot(sol_df); colors = colors)
    end
end

function create_beta_plot(beta_df)
    return data(beta_df) *
           mapping(:time => "Time (days)", :β) *
           visual(Lines; linewidth = 4, color = "green")
end

function draw_beta_plot(beta_plot)
    return draw(beta_plot)
end

function draw_beta_plot(beta_df::DataFrame)
    return draw_beta_plot(create_beta_plot(beta_df))
end

function draw_combined_sir_beta_plot(sir_plot, beta_plot)
    combined = Figure()

    axsir = Axis(combined[2:4, 1]; xlabel = "Time (days)", ylabel = "Number of individuals")
    axbeta = Axis(
        combined[1, 1]; ylabel = "β", xticksvisible = false, xticklabelsvisible = false
    )

    draw!(axsir, sir_plot)
    draw!(axbeta, beta_plot)

    return combined
end

function create_sir_quantiles_plot(
    sim_quantiles = sim_quantiles; lower, upper, quantiles, δt = δt
)
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
        color = colors[1],
        linewidth = 2,
        label = "S",
    )
    lines!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[med_index, :, 2];
        color = colors[2],
        linewidth = 2,
        label = "I",
    )
    lines!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[med_index, :, 3];
        color = colors[3],
        linewidth = 2,
        label = "R",
    )
    lines!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[med_index, :, 4];
        color = colors[4],
        linewidth = 2,
        label = "N",
    )

    # User-specified quantiles
    band!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[lower_index, :, 1],
        sim_quantiles[upper_index, :, 1];
        color = (colors[1], 0.5),
    )
    band!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[lower_index, :, 2],
        sim_quantiles[upper_index, :, 2];
        color = (colors[2], 0.5),
    )
    band!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[lower_index, :, 3],
        sim_quantiles[upper_index, :, 3];
        color = (colors[3], 0.5),
    )
    band!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[lower_index, :, 4],
        sim_quantiles[upper_index, :, 4];
        color = (colors[4], 0.5),
    )

    Legend(fig[1, 2], ax, "State")

    return fig
end

function create_beta_quantiles_plot(
    sim_quantiles = sim_quantiles; lower, upper, quantiles, δt = δt
)
    fig = Figure()
    ax = Axis(fig[1, 1])

    med_index = findfirst(isequal(0.5), quantiles)
    lower_index = findfirst(isequal(lower), quantiles)
    upper_index = findfirst(isequal(upper), quantiles)

    # Medians
    lines!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[med_index, :, 5];
        color = colors[1],
        linewidth = 2,
        label = "β",
    )

    # User-specified quantiles
    band!(
        ax,
        tlower:δt:tmax,
        sim_quantiles[lower_index, :, 5],
        sim_quantiles[upper_index, :, 5];
        color = ("green", 0.5),
    )

    Legend(fig[1, 2], ax, "Parameter")

    return fig
end