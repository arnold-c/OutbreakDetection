"""
This file contains functions to help visually inspect ensemble simulations
"""

function incidence_prevalence_plot(
    incidencearr,
    ensemblearr,
    thresholdsarr,
    timeparams;
    colormap = [N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR],
    threshold = 5,
    kwargs...,
)
    @unpack trange = timeparams
    times = collect(trange) ./ 365
    kwargs_dict = Dict(kwargs)

    period_sum_arr = zeros(Int64, length(times), 2)
    for (lower, upper, periodsum, outbreakstatus) in
        eachrow(thresholdsarr[1])
        period_sum_arr[lower:upper, 1] .= periodsum
        period_sum_arr[lower:upper, 2] .= outbreakstatus
    end

    fig = Figure()
    ax_prev = Axis(fig[1, 1]; ylabel = "Prevalence")
    ax_inc = Axis(fig[2, 1]; ylabel = "Incidence")
    ax_periodsum = Axis(
        fig[3, 1]; xlabel = "Time (years)", ylabel = "Period Sum"
    )

    linkxaxes!(ax_prev, ax_inc, ax_periodsum)

    lines!(
        ax_prev,
        times,
        ensemblearr[:, 2, 1];
        color = period_sum_arr[:, 2],
        colormap = colormap,
    )
    lines!(
        ax_inc,
        times,
        incidencearr[:, 1, 1];
        color = period_sum_arr[:, 2],
        colormap = colormap,
    )
    hlines!(
        ax_inc,
        threshold;
        color = :black,
        linestyle = :dash,
        linewidth = 2,
    )
    barplot!(
        ax_periodsum,
        times,
        period_sum_arr[:, 1];
        color = period_sum_arr[:, 2],
        colormap = colormap,
    )

    map(hidexdecorations!, [ax_prev, ax_inc])

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [ax_prev, ax_inc, ax_periodsum],
        )
    end

    if haskey(kwargs_dict, :ylims_periodsum)
        ylims!(ax_periodsum, kwargs_dict[:ylims_periodsum])
    end

    if haskey(kwargs_dict, :ylims_inc)
        ylims!(ax_inc, kwargs_dict[:ylims_inc])
    end

    Legend(
        fig[:, 2],
        [PolyElement(; color = col) for col in colormap],
        ["Not Outbreak", "Outbreak"],
        "True\nOutbreak Status",
    )

    return fig
end

function incidence_testing_plot(
    incarr,
    noisearr,
    testingarr,
    test_movingvg_arr,
    detection_specification,
    timeparams;
    sim = 1,
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR, N_ALERTS_COLOR
    ],
    plottitle = "",
    kwargs...,
)
    times = collect(timeparams.trange) ./ 365
    kwargs_dict = Dict(kwargs)

    inc_test_fig = Figure()
    inc_test_ax1 = Axis(inc_test_fig[1, 1]; ylabel = "Incidence")
    inc_test_ax2 = Axis(inc_test_fig[2, 1]; ylabel = "Incidence + Noise")
    inc_test_ax3 = Axis(inc_test_fig[3, 1]; ylabel = "Test Positive")
    inc_test_ax4 = Axis(
        inc_test_fig[4, 1];
        xlabel = "Time (years)",
        ylabel = "7d Avg Test Positive",
    )

    lines!(
        inc_test_ax1, times, incarr[:, 1, sim];
        color = incarr[:, 3, sim],
        colormap = outbreakcolormap,
    )
    lines!(
        inc_test_ax2, times, incarr[:, 1, sim] .+ noisearr[:, sim];
        color = incarr[:, 3, sim],
        colormap = outbreakcolormap,
    )
    lines!(
        inc_test_ax3, times, testingarr[:, 5, sim];
        color = testingarr[:, 6, sim],
        colormap = alertcolormap,
    )
    lines!(
        inc_test_ax4, times, test_movingvg_arr[:, sim];
        color = testingarr[:, 6, sim],
        colormap = alertcolormap,
    )

    linkxaxes!(inc_test_ax1, inc_test_ax2, inc_test_ax3, inc_test_ax4)

    map(hidexdecorations!, [inc_test_ax1, inc_test_ax3])

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [inc_test_ax1, inc_test_ax2, inc_test_ax3, inc_test_ax4],
        )
    end

    if haskey(kwargs_dict, :ylims)
        map(
            ax -> ylims!(ax, kwargs_dict[:ylims]),
            [inc_test_ax1, inc_test_ax2, inc_test_ax3, inc_test_ax4],
        )
    end

    map(
        ax -> hlines!(
            ax,
            5;
            color = :black,
            linestyle = :dash,
            linewidth = 2,
        ),
        [inc_test_ax1, inc_test_ax2],
    )

    map(
        ax -> hlines!(
            ax,
            detection_specification.alert_threshold;
            color = :black,
            linestyle = :dash,
            linewidth = 2,
        ),
        [inc_test_ax3, inc_test_ax4],
    )

    Label(
        inc_test_fig[0, :, Top()],
        plottitle,
    )

    Legend(
        inc_test_fig[1:2, 2],
        [PolyElement(; color = col) for col in outbreakcolormap],
        ["Not Outbreak", "Outbreak"],
        "True\nOutbreak Status",
    )

    Legend(
        inc_test_fig[3:4, 2],
        [PolyElement(; color = col) for col in alertcolormap],
        ["Not Outbreak", "Outbreak"],
        "Alert Status",
    )

    rowsize!(inc_test_fig.layout, 0, 5)
    colsize!(inc_test_fig.layout, 1, Relative(0.92))

    return inc_test_fig
end

function testing_plot(
    testingarr, timeparams; plottitle = "", sim1_num = 1, sim2_num = 25
)
    times = collect(timeparams.trange) ./ 365

    testing_fig = Figure()
    testing_grid = testing_fig[1, 1] = GridLayout()
    sim1_ax = Axis(
        testing_grid[1, 1];
        title = "Simulation $sim1_num",
        xlabel = "Time (years)",
        ylabel = "Number tested",
    )
    sim2_ax = Axis(
        testing_grid[2, 1];
        title = "Simulation $sim2_num",
        xlabel = "Time (years)",
        ylabel = "Number tested",
    )

    for (ind, label, col) in zip(
        (1, 2, 5),
        ("Infectious", "Noise", "Total Positive"),
        (:red, :blue, :black),
    )
        lines!(
            sim1_ax, times, testingarr[:, ind, sim1_num]; color = col,
            label = label,
        )
        lines!(
            sim2_ax, times, testingarr[:, ind, sim2_num]; color = col,
            label = label,
        )
    end

    linkxaxes!(sim1_ax, sim2_ax)

    hidexdecorations!(sim1_ax)

    Legend(
        testing_fig[:, 2],
        sim1_ax,
        "Type of Individual";
        orientation = :vertical,
    )

    Label(
        testing_fig[0, :, Top()],
        plottitle,
    )

    rowsize!(testing_fig.layout, 0, 5)
    colsize!(testing_fig.layout, 1, Relative(0.92))

    return testing_fig
end

function ensemble_outbreak_distribution_plot(testarr, infecarr; plottitle = "")
    outbreak_dist_fig = Figure()
    outbreak_dist_ax = Axis(
        outbreak_dist_fig[1, 1];
        xlabel = "Proportion of Time Series with Outbreak",
    )

    hist!(
        outbreak_dist_ax,
        vec(sum(@view(infecarr[:, 3, :]); dims = 1)) ./ size(infecarr, 1);
        bins = 0.0:0.01:0.7,
        color = (:blue, 0.5),
        strokecolor = :black,
        strokewidth = 1,
        normalization = :pdf,
        label = "True Outbreaks",
    )

    hist!(
        outbreak_dist_ax,
        vec(sum(@view(testarr[:, 6, :]); dims = 1)) ./ size(testarr, 1);
        bins = 0.0:0.01:0.7,
        color = (:red, 0.5),
        strokecolor = :black,
        strokewidth = 1,
        normalization = :pdf,
        label = "Tested Outbreaks",
    )

    Legend(outbreak_dist_fig[1, 2], outbreak_dist_ax, "Outbreak Proportion")

    Label(
        outbreak_dist_fig[0, :, Top()],
        plottitle,
    )

    rowsize!(outbreak_dist_fig.layout, 0, 5)
    colsize!(outbreak_dist_fig.layout, 1, Relative(0.92))

    return outbreak_dist_fig
end
