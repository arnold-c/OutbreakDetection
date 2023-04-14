function create_sir_plot(sol_df; labels = ["S", "I", "R", "N"], annual = annual)
    time_function(t) = annual ==true ? t / 365.0 : t
    if annual == true
        time_label = "Time (years)"
    else
        time_label = "Time (days)"
    end

    return data(sol_df) *
           mapping(
                :time => time_function => time_label, :Number;
               color = :State => sorter(labels...),
           ) *
           visual(Lines; linewidth = 4)
end

function draw_sir_plot(
    sir_plot; colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"], xlims = xlims, ylims = ylims,
)
    return draw(sir_plot; palettes = (; color = colors), axis = (; limits = (xlims, ylims)),)
end

function draw_sir_plot(
    sol_df::DataFrame;
    colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"],
    labels = ["S", "I", "R", "N"],
    annual = false,
    xlims = nothing,
    ylims = nothing,
)
    return draw_sir_plot(
        create_sir_plot(sol_df; labels = labels, annual = annual);
        colors = colors,
        xlims = xlims,
        ylims = ylims
    )
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

    axsir = Axis(
        combined[2:4, 1];
        xlabel = "Time (days)",
        ylabel = "Number of individuals",
    )
    axbeta = Axis(
        combined[1, 1]; ylabel = "β", xticksvisible = false,
        xticklabelsvisible = false,
    )

    draw!(axsir, sir_plot)
    draw!(axbeta, beta_plot)

    return combined
end

function sir_quantiles_array_base_plot(
    sim_quantiles, lower_index, med_index, upper_index, δt, colors, labels,
    annual; xlims, ylims, caption,
)
    times = tlower:δt:tmax
    xlab = "Time (days)"
    if annual == true
        times = times ./ 365
        xlab = "Time (years)"
    end

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = xlab, ylabel = "Number")

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
        eachindex(labels),
    )

    map(
        state -> band!(
            ax,
            times,
            sim_quantiles[lower_index, :, state],
            sim_quantiles[upper_index, :, state];
            color = (colors[state], 0.5),
        ),
        eachindex(labels),
    )

    if xlims != false
        xlims!(ax, xlims)
    end

    if ylims != false
        ylims!(ax, ylims)
    end

    Legend(fig[1, 2], ax, "State")

    if caption != false
        Label(fig[2, :, Bottom()], caption)
        rowsize!(fig.layout, 1, Relative(0.98))
    end

    return fig
end

function create_sir_quantiles_plot(
    sim_quantiles = sim_quantiles;
    δt = δt,
    colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"],
    labels = ["S", "I", "R", "N"],
    annual = false,
    xlims = false,
    ylims = false,
    caption = false,
)
    med_index = 2
    lower_index = 1
    upper_index = 3

    return sir_quantiles_array_base_plot(
        sim_quantiles, lower_index, med_index, upper_index, δt,
        colors,
        labels,
        annual;
        xlims = xlims,
        ylims = ylims,
        caption = caption,
    )
end

function create_sir_quantiles_plot(
    sim_quantiles, lower, upper, quantiles;
    δt = δt,
    colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"],
    labels = ["S", "I", "R", "N"],
    annual = false,
)
    med_index = findfirst(isequal.(0.5), quantiles)
    lower_index = findfirst(isequal.(lower), quantiles)
    upper_index = findfirst(isequal.(upper), quantiles)

    return sir_quantiles_array_base_plot(
        sim_quantiles, lower_index, med_index, upper_index, δt,
        colors, labels, annual,
    )
end

function create_sir_quantiles_plot(
    summ::EnsembleSummary; annual = false,
    colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"],
    labels = ["S", "I", "R", "N"],
)
    times = tlower:δt:tmax
    xlab = "Time (days)"
    if annual == true
        times = times ./ 365
        xlab = "Time (years)"
    end

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = xlab, ylabel = "Number")

    med_df = DataFrame(summ.med)[:, 2:end]
    lower = DataFrame(summ.qlow)[:, 2:end]
    upper = DataFrame(summ.qhigh)[:, 2:end]

    map(
        x -> transform!(x, Cols(:) => (+) => :N),
        [med_df, lower, upper],
    )

    map(
        state -> lines!(
            ax,
            times,
            med_df[:, state];
            color = colors[state],
            linewidth = 2,
            label = labels[state],
        ),
        eachindex(labels),
    )

    map(
        state -> band!(
            ax,
            times,
            lower[:, state],
            upper[:, state];
            color = (colors[state], 0.5),
        ),
        eachindex(labels),
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