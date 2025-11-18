export plot_schematic

function plot_schematic(
        inc_vec,
        outbreakstatus_vec,
        outbreak_bounds,
        outbreak_specification,
        noise_vec,
        testpositive_vec,
        alertstatus_vec,
        alert_bounds,
        alertthreshold;
        time_p = SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 10, tstep = 1.0),
        outbreakcolormap = [
            "#5E5C6C",
            "#CE6F58",
        ],
        alertcolormap = [
            "#F4A157",
            "#2B3465",
        ],
        shade_alert_outbreak_overlap = false,
        measlesalpha = 0.5,
        testalpha = 0.5,
        kwargs...,
    )
    kwargs_dict = Dict(kwargs)

    times = collect(time_p.trange)

    if haskey(kwargs_dict, :xlims)
        lower = maximum([1, kwargs_dict[:xlims][1] * 365])
        upper = minimum([Int64(time_p.tlength), kwargs_dict[:xlims][2] * 365])
        times = times[lower:upper]
        inc_vec = inc_vec[lower:upper]
        outbreakstatus_vec = outbreakstatus_vec[lower:upper]
        noise_vec = noise_vec[lower:upper]
        testpositive_vec = testpositive_vec[lower:upper]
        alertstatus_vec = alertstatus_vec[lower:upper]

        outbreak_bounds = outbreak_bounds[
            (outbreak_bounds[:, 1] .>= lower) .& (outbreak_bounds[:, 2] .<= upper),
            :,
        ]

        alert_bounds = alert_bounds[
            (alert_bounds[:, 1] .>= lower) .& (alert_bounds[:, 2] .<= upper), :,
        ]
    end
    outbreak_bounds_vec = vec(outbreak_bounds[:, 1:2])
    alert_bounds_vec = vec(alert_bounds[:, 1:2])

    fig = Figure()
    incga = fig[1, 1] = GridLayout()
    noisega = fig[2, 1] = GridLayout()
    testga = fig[3, 1] = GridLayout()
    incax = Axis(incga[1, 1]; ylabel = "Measles Incidence")
    noiseax = Axis(noisega[1, 1]; ylabel = "Noise Incidence")
    testax = Axis(testga[1, 1]; xlabel = "Time", ylabel = "Test Positives")

    if shade_alert_outbreak_overlap
        if !isempty(outbreak_bounds_vec)
            map(
                ax -> vspan!(
                    ax,
                    outbreak_bounds[:, 1],
                    outbreak_bounds[:, 2];
                    color = (outbreakcolormap[2], measlesalpha),
                ),
                [incax, noiseax, testax],
            )
        end

        if !isempty(alert_bounds)
            vspan!(
                testax,
                alert_bounds[:, 1],
                alert_bounds[:, 2];
                color = (alertcolormap[2], testalpha),
            )
        end
    end

    lines!(
        incax,
        times,
        inc_vec;
        color = outbreakstatus_vec,
        colormap = outbreakcolormap,
        linewidth = 3,
    )

    hlines!(
        incax,
        outbreak_specification.outbreak_threshold;
        color = :black,
        linewidth = 2,
        linestyle = :dash,
    )

    lines!(
        noiseax,
        times,
        noise_vec;
        color = outbreakcolormap[1],
        linewidth = 3,
    )

    lines!(
        testax,
        times,
        testpositive_vec;
        color = alertstatus_vec,
        colormap = alertcolormap,
        linewidth = 3,
    )

    hlines!(
        testax,
        alertthreshold;
        color = :black,
        linewidth = 2,
        linestyle = :dash,
    )

    text!(
        testax,
        times[1] + 1,
        alertthreshold + 0.5;
        text = "T = $alertthreshold",
        justification = :left,
    )

    for (label, layout) in
        zip(["a", "b", "c"], [incga[1, 1], noisega[1, 1], testga[1, 1]])
        Label(
            layout[1, 1, TopLeft()], label;
            fontsize = 30,
            font = :bold,
            padding = (0, 0, 20, 0),
            halign = :right
        )
    end

    map(
        ax -> hidexdecorations!(ax),
        [
            noiseax,
            incax,
            testax,
        ],
    )

    linkxaxes!(noiseax, incax, testax)

    return fig
end
