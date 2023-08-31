using DrWatson
@quickactivate "OutbreakDetection"

using GLMakie
using AlgebraOfGraphics
using ColorSchemes
using UnPack
using DataFrames

set_aog_theme!()
# Set depending on size of screen
update_theme!(; resolution = (2200, 1300))
# update_theme!(; resolution = (850, 600))
GLMakie.activate!(; float = true)

seircolors = ["dodgerblue4", "green", "firebrick3", "chocolate2", "purple"]
seir_state_labels = ["S", "E", "I", "R", "N"]

function create_sir_plot(sol_df; labels = ["S", "I", "R", "N"], annual = annual)
    time_function(t) = annual == true ? t / 365.0 : t
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
    sir_plot; colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"],
    xlims = xlims, ylims = ylims,
)
    return draw(
        sir_plot;
        palettes = (; color = colors),
        axis = (; limits = (xlims, ylims)),
    )
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
        ylims = ylims,
    )
end

function sir_quantiles_array_base_plot(
    sim_quantiles, lower_index, med_index, upper_index, timeparams, colors,
    labels,
    annual; xlims, ylims, caption,
)
    times = timeparams.trange
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
    timeparams,
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
        sim_quantiles, lower_index, med_index, upper_index, timeparams,
        colors,
        labels,
        annual;
        xlims = xlims,
        ylims = ylims,
        caption = caption,
    )
end

outbreakcols = [ColorSchemes.magma[i] for i in (200, 20)]

function detect_outbreak_plot(
    incidencearr, ensemblearr, timeparams; colormap = outbreakcols, kwargs...
)
    @unpack trange = timeparams
    times = collect(trange) ./ 365
    kwargs_dict = Dict(kwargs)

    fig = Figure()
    ax_prev = Axis(fig[1, 1]; ylabel = "Prevalence")
    ax_inc = Axis(fig[2, 1]; ylabel = "Incidence")
    ax_periodsum = Axis(
        fig[3, 1]; xlabel = "Time (years)", ylabel = "Period Sum"
    )

    linkxaxes!(ax_prev, ax_inc, ax_periodsum)

    lines!(ax_prev, times, ensemblearr[2, :, 1])
    lines!(ax_inc, times, incidencearr[:, 1, 1])
    barplot!(
        ax_periodsum,
        times,
        incidencearr[:, 3, 1];
        color = incidencearr[:, 4, 1],
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

    axislegend(
        ax_periodsum,
        [PolyElement(; color = col) for col in colormap],
        ["Not Outbreak", "Outbreak"];
        bgcolor = :white,
        framecolor = :white,
        framevisible = true,
        padding = (10.0f0, 10.0f0, 8.0f0, 8.0f0),
    )

    return fig
end

function visualize_ensemble_noise(ensemble_noise_arr, timeparams)
    times = collect(timeparams.trange) ./ 365
    noise_fig = Figure()
    noise_ax = Axis(
        noise_fig[1, 1]; xlabel = "Time (years)", ylabel = "Noise Incidence"
    )

    for sim in axes(ensemble_noise_arr, 3)
        lines!(
            noise_ax,
            times,
            ensemble_noise_arr[:, 1, sim];
            color = (:red, 0.1)
        )
    end

    return noise_fig
end

function incidence_testing_plot(
    incarr,
    testingarr,
    timeparams,
    detectthreshold;
    sim = 1,
    colormap = outbreakcols,
    kwargs...,
)
    times = collect(timeparams.trange) ./ 365
    kwargs_dict = Dict(kwargs)

    inc_test_fig = Figure()
    inc_test_ax1 = Axis(inc_test_fig[1, 1]; ylabel = "Incidence")
    inc_test_ax2 = Axis(inc_test_fig[2, 1]; ylabel = "Test Positive")
    inc_test_ax3 = Axis(
        inc_test_fig[3, 1];
        xlabel = "Time (years)",
        ylabel = "7d Avg Test Positive",
    )

    lines!(
        inc_test_ax1, times, incarr[:, 1, sim];
        color = incarr[:, 4, sim],
        colormap = colormap,
    )
    lines!(
        inc_test_ax2, times, testingarr[:, 3, sim];
        color = testingarr[:, 7, sim],
        colormap = colormap,
    )
    lines!(
        inc_test_ax3, times, testingarr[:, 6, sim];
        color = testingarr[:, 7, sim],
        colormap = colormap,
    )

    linkxaxes!(inc_test_ax1, inc_test_ax2, inc_test_ax3)

    map(hidexdecorations!, [inc_test_ax1, inc_test_ax2])

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [inc_test_ax1, inc_test_ax2, inc_test_ax3],
        )
    end

    if haskey(kwargs_dict, :ylims)
        map(
            ax -> ylims!(ax, kwargs_dict[:ylims]),
            [inc_test_ax1, inc_test_ax2, inc_test_ax3],
        )
    end

    hlines!(
        inc_test_ax1, 5;
        color = :black,
        linestyle = :dash,
        linewidth = 2,
    )

    map(
        ax -> hlines!(
            ax,
            detectthreshold;
            color = :black,
            linestyle = :dash,
            linewidth = 2,
        ),
        [inc_test_ax2, inc_test_ax3],
    )

    Legend(
        inc_test_fig[1, 2],
        [PolyElement(; color = col) for col in outbreakcols],
        ["Not Outbreak", "Outbreak"],
        "True\nOutbreak Status",
    )

    Legend(
        inc_test_fig[2:3, 2],
        [PolyElement(; color = col) for col in outbreakcols],
        ["Not Outbreak", "Outbreak"],
        "Detected\nOutbreak Status",
    )

    return inc_test_fig
end
