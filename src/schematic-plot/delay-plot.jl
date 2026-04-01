export trim_schematic_window,
    plot_test_delay_panel!

function trim_schematic_window(
        values,
        outbreak_bounds,
        alert_bounds,
        time_parameters;
        xlims,
    )
    lower = maximum([1, xlims[1] * 365])
    upper = minimum([Int64(time_parameters.tlength), xlims[2] * 365])
    time_window = lower:upper

    trimmed_outbreak_bounds = outbreak_bounds[
        (outbreak_bounds[:, 1] .>= lower) .& (outbreak_bounds[:, 2] .<= upper), 1:2,
    ]
    trimmed_alert_bounds = alert_bounds[
        (alert_bounds[:, 1] .>= lower) .& (alert_bounds[:, 2] .<= upper), :,
    ]

    return (
        time_window = time_window,
        values = values[time_window],
        outbreak_bounds = trimmed_outbreak_bounds,
        alert_bounds = trimmed_alert_bounds,
    )
end


function plot_test_delay_panel!(
        ax,
        label,
        panel_data,
        alert_status,
        time_parameters,
        alert_threshold;
        alertcolormap,
        outbreakcolormap,
        measlesalpha,
        testalpha,
        linewidth,
        ylabelsize,
        thresholdfontsize,
        xticklabelsize,
        yticklabelsize,
        show_xlabel = false,
    )
    times = collect(time_parameters.trange)[panel_data.time_window]
    year_ticks = ceil(Int, first(times) / 365):floor(Int, last(times) / 365)

    if !isempty(panel_data.outbreak_bounds)
        vspan!(
            ax,
            panel_data.outbreak_bounds[:, 1],
            panel_data.outbreak_bounds[:, 2];
            color = (outbreakcolormap[2], measlesalpha),
        )
    end

    if !isempty(panel_data.alert_bounds)
        vspan!(
            ax,
            panel_data.alert_bounds[:, 1],
            panel_data.alert_bounds[:, 2];
            color = (alertcolormap[2], testalpha),
        )
    end

    lines!(
        ax,
        times,
        panel_data.values;
        color = alert_status[panel_data.time_window],
        colormap = alertcolormap,
        linewidth = linewidth,
    )

    hlines!(
        ax,
        alert_threshold;
        color = :black,
        linewidth = linewidth,
        linestyle = :dash,
    )

    text!(
        ax,
        times[1] + 1,
        alert_threshold + 0.5;
        text = "T = $alert_threshold",
        justification = :left,
        fontsize = thresholdfontsize,
    )

    ax.title = label
    ax.ylabel = "Test Positives"
    ax.ylabelsize = ylabelsize
    ax.xticklabelsize = xticklabelsize
    ax.yticklabelsize = yticklabelsize
    ax.xticks = (collect(year_ticks .* 365), string.(collect(year_ticks)))

    return if show_xlabel
        ax.xlabel = "Time (years)"
    else
        hidexdecorations!(ax; grid = false, label = true)
    end
end
