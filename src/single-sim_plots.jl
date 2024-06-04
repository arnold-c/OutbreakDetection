seircolors = ["dodgerblue4", "green", "firebrick3", "chocolate2", "purple"]
seir_state_labels = ["S", "E", "I", "R", "N"]

function single_seir_plot(
    sol_df;
    state_colors = seircolors,
    state_labels = seir_state_labels,
    annual = true,
    kwargs...,
)
    times = time_function.(unique(sol_df.time); annual = annual)
    time_label = annual_label(; annual)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = time_label, ylabel = "Number")

    for (i, label) in pairs(state_labels)
        lines!(
            ax,
            times,
            subset(sol_df, :State => s -> s .== label)[!, "Number"];
            label = label,
            color = state_colors[i],
            linewidth = 4,
        )
    end

    lims_check(ax, kwargs)

    Legend(fig[1, 2], ax, "State")

    return fig
end

function single_seir_statespace_plot(
    seir_array;
    colormap = :magma,
    annual = true,
    kwargs...,
)
    times = time_function.(axes(seir_array, 1); annual = annual)
    time_label = annual_label(; annual)

    fig = Figure()

    ax = Axis(fig[1, 1]; xlabel = "I", ylabel = "S")

    lines!(
        ax,
        seir_array[:, 3],
        seir_array[:, 1];
        color = axes(seir_array, 1),
        colormap = colormap,
    )

    lims_check(ax, kwargs)

    Colorbar(
        fig[1, 2];
        colormap = colormap,
        limits = (minimum(times), maximum(times)),
        label = time_label,
    )

    return fig
end

function time_function(t; annual = true)
    return annual == true ? t / 365.0 : t
end

function annual_label(; annual = true)
    return annual == true ? "Time (years)" : "Time (days)"
end

function lims_check(ax, kwargs)
    if haskey(kwargs, :xlims)
        xlims!(ax, kwargs[:xlims])
    end

    if haskey(kwargs, :ylims)
        ylims!(ax, kwargs[:ylims])
    end
end
