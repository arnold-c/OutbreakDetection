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
    sim_quantiles = sim_quantiles; lower = 0.1, upper = 0.9, quantiles = quantiles, δt = δt,
    colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"], labels = ["S", "I", "R", "N"], 
    annual = false,
)
    times = tlower:δt:tmax
    xlab = "Time (days)"
    if annual == true
        times = times ./ 365
        xlab = "Time (years)"
    end

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = xlab, ylabel = "Number")

    med_index = findfirst(isequal(0.5), quantiles)
    lower_index = findfirst(isequal(lower), quantiles)
    upper_index = findfirst(isequal(upper), quantiles)

    # Medians
    map(
        state -> lines!(
            ax,
            times,
            sim_quantiles[med_index, :, state];
            color = colors[state],
            linewidth = 2,
            label = labels[state], 
        ),
        eachindex(labels)
    )

    map(
        state -> band!(
            ax,
            times,
            sim_quantiles[lower_index, :, state],
            sim_quantiles[upper_index, :, state];
            color = (colors[state], 0.5)
        ),
        eachindex(labels)  
    )

    Legend(fig[1, 2], ax, "State")

    return fig
end

function create_sir_quantiles_plot(
    summ::EnsembleSummary; annual = false,
    colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"],
    labels = ["S", "I", "R", "N"]
)

    times = tlower:δt:tmax
    xlab = "Time (days)"
    if annual == true
        times = times ./ 365
        xlab = "Time (years)"
    end

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = xlab, ylabel = "Number")

    med_df = DataFrame(summ.med)[:, 2:end]
    lower = DataFrame(summ.qlow)[:, 2:end]
    upper = DataFrame(summ.qhigh)[:, 2:end]

    map(
        x -> transform!(x, Cols(:) => (+) => :N),
        [med_df, lower, upper]
    )

    map(
        state -> lines!(
            ax,
            times,
            med_df[:, state],
            color = colors[state],
            linewidth = 2,
            label = labels[state]
        ),
        eachindex(labels)
    )

    map(
        state -> band!(
            ax,
            times,
            lower[:, state],
            upper[:, state],
            color = (colors[state], 0.5)
        ),
        eachindex(labels)
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

    times = tlower:δt:tmax
    if annual == true
        times = times ./ 365
    end

    # Medians
    lines!(
        ax,
        times,
        sim_quantiles[med_index, :, 5];
        color = colors[1],
        linewidth = 2,
        label = "β",
    )

    # User-specified quantiles
    band!(
        ax,
        times,
        sim_quantiles[lower_index, :, 5],
        sim_quantiles[upper_index, :, 5];
        color = ("green", 0.5),
    )

    Legend(fig[1, 2], ax, "Parameter")

    return fig
end