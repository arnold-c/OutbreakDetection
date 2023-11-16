# module PlottingFunctions
#
# export seircolors, seir_state_labels, create_sir_plot, draw_sir_plot,
#     sir_quantiles_array_base_plot, create_sir_quantiles_plot, outbreakcols,
#     detect_outbreak_plot, visualize_ensemble_noise, incidence_testing_plot,
#     testing_plot, ensemble_outbreak_distribution_plot, ensemble_OTChars_plot

using GLMakie
using AlgebraOfGraphics
using ColorSchemes
using UnPack
using DataFrames
using FLoops
using LaTeXStrings
using NaNMath: NaNMath

seircolors = ["dodgerblue4", "green", "firebrick3", "chocolate2", "purple"]
seir_state_labels = ["S", "E", "I", "R", "N"]

ACCURACY_COLOR = "#004643"
DAILY_SENSITIVITY_COLOR = "#2A3965"
DAILY_SPECIFICITY_COLOR = "#C31D60"
DAILY_PPV_COLOR = "#22D37D"
DAILY_NPV_COLOR = "#6f366bff"
DETECTION_DELAY_COLOR = "#AE560A"
N_ALERTS_PER_OUTBREAK_COLOR = "#86B1A3"
N_FALSE_ALERTS_COLOR = "#D06778"
N_ALERTS_COLOR = "#00857E"
N_OUTBREAKS_COLOR = "#F4A157"
N_MISSED_OUTBREAKS_COLOR = "#5E5C6C"
PERC_OUTBREAKS_DETECTED_COLOR = "#F0780F"
PERC_OUTBREAKS_MISSED_COLOR = "#3A3842"
PERC_ALERTS_CORRECT_COLOR = "#004643"
PERC_ALERTS_FALSE_COLOR = "#852938"

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

function bifurcation_plot(
    x_vector,
    annual_summary;
    years,
    xlabel = "Birth rate (per 1_000, per annum)",
    ylabel = "Max. I",
)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel)

    for year in eachindex(years)
        scatter!(
            ax,
            x_vector,
            annual_summary[year, 2, :];
            markersize = 4,
            color = :black,
        )
    end

    return fig
end

function bifurcation_heatmap(
    birth_rate_vec,
    beta_force_vec,
    cycle_summary,
)
    fig, ax, hm = heatmap(
        birth_rate_vec,
        beta_force_vec,
        cycle_summary,
    )

    Colorbar(
        fig[:, end + 1],
        hm;
        label = "Periodicity",
    )

    ax.xlabel = "Birth rate (per 1_000, per annum)"
    ax.ylabel = "beta_force (seasonality)"

    return fig
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

function visualize_ensemble_noise(
    ensemble_inc_arr, ensemble_noise_spec, timeparams
)
    ensemble_noise_arr = create_poisson_noise_arr(
        ensemble_inc_arr, ensemble_noise_spec
    )

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
    noisearr,
    testingarr,
    timeparams,
    alertthreshold;
    sim = 1,
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR, N_ALERTS_COLOR
    ],
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
        inc_test_ax2, times, incarr[:, 1, sim] .+ noisearr[:, 1, sim];
        color = incarr[:, 3, sim],
        colormap = outbreakcolormap,
    )
    lines!(
        inc_test_ax3, times, testingarr[:, 3, sim];
        color = testingarr[:, 7, sim],
        colormap = alertcolormap,
    )
    lines!(
        inc_test_ax4, times, testingarr[:, 6, sim];
        color = testingarr[:, 7, sim],
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
            alertthreshold;
            color = :black,
            linestyle = :dash,
            linewidth = 2,
        ),
        [inc_test_ax3, inc_test_ax4],
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

    return inc_test_fig
end

function testing_plot(testingarr, timeparams)
    times = collect(timeparams.trange) ./ 365

    testing_fig = Figure()
    testing_grid = testing_fig[1, 1] = GridLayout()
    sim1_ax = Axis(
        testing_grid[1, 1];
        title = "Simulation 1",
        xlabel = "Time (years)",
        ylabel = "Number tested",
    )
    sim2_ax = Axis(
        testing_grid[2, 1];
        title = "Simulation 2",
        xlabel = "Time (years)",
        ylabel = "Number tested",
    )
    sim3_ax = Axis(
        testing_grid[1, 2];
        title = "Simulation 3",
        xlabel = "Time (years)",
        ylabel = "Number tested",
    )
    sim4_ax = Axis(
        testing_grid[2, 2];
        title = "Simulation 4",
        xlabel = "Time (years)",
        ylabel = "Number tested",
    )

    for (ind, label, col) in zip(
        (1, 2, 5),
        ("Infectious", "Noise", "Total Positive"),
        (:red, :blue, :black),
    )
        lines!(
            sim1_ax, times, testingarr[:, ind, 1]; color = col, label = label
        )
        lines!(
            sim2_ax, times, testingarr[:, ind, 2]; color = col, label = label
        )
        lines!(
            sim3_ax, times, testingarr[:, ind, 3]; color = col, label = label
        )
        lines!(
            sim4_ax, times, testingarr[:, ind, 4]; color = col, label = label
        )
    end

    linkxaxes!(sim1_ax, sim2_ax)
    linkxaxes!(sim3_ax, sim4_ax)

    linkyaxes!(sim1_ax, sim3_ax)
    linkyaxes!(sim2_ax, sim4_ax)

    map(hidexdecorations!, [sim1_ax, sim3_ax])
    map(hideydecorations!, [sim3_ax, sim4_ax])

    Legend(
        testing_fig[2, :],
        sim1_ax,
        "Type of Individual";
        orientation = :horizontal,
    )

    return testing_fig
end

function ensemble_outbreak_distribution_plot(testarr, infecarr)
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
        vec(sum(@view(testarr[:, 7, :]); dims = 1)) ./ size(testarr, 1);
        bins = 0.0:0.01:0.7,
        color = (:red, 0.5),
        strokecolor = :black,
        strokewidth = 1,
        normalization = :pdf,
        label = "Tested Outbreaks",
    )

    Legend(outbreak_dist_fig[1, 2], outbreak_dist_ax, "Outbreak Proportion")

    return outbreak_dist_fig
end

function ensemble_OTChars_plot(
    OTChars,
    testspec,
    detectspec,
    plottingchars;
    columnfacetchar = :alert_threshold,
    columnfacetchar_label = "Alert Threshold",
    xlabel = "Alert Characteristic Value",
    ylabel = "Density",
    binwidth = 10.0,
    legend = true,
    legendlabel = "# Outbreaks",
    meanlines = true,
    meanlabels = true,
    normalization = :none,
    kwargs...,
)
    otchars_vec = Vector{NamedTuple}(undef, 1)
    otchars_vec[1] = (;
        OT_chars = OTChars,
        ind_test_spec = testspec,
        outbreak_detect_spec = detectspec,
    )
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    @pack! kwargs_dict = columnfacetchar_label,
    binwidth, xlabel, ylabel,
    legendlabel

    fig = Figure()

    construct_OTchars_facets!(
        fig,
        otchars_vec,
        plottingchars,
        [1],
        [1],
        columnfacetchar,
        kwargs_dict;
        meanlines = meanlines,
        meanlabels = meanlabels,
        normalization = normalization,
    )

    if legend
        Legend(
            fig[1, 2],
            [
                PolyElement(; color = col) for
                col in map(chartuple -> chartuple.color, plottingchars)
            ],
            [chartuple.label for chartuple in plottingchars];
            label = legendlabel,
        )
    end

    return fig
end

function ensemble_outbreak_detect_diff_plot(OT_chars; binwidth = 1)
    fig = Figure()
    ax = Axis(
        fig[1, 1]; xlabel = "Difference Between Actual - Detected Outbreaks"
    )

    difference = OT_chars.noutbreaks .- OT_chars.nalerts

    bins = minimum(difference):binwidth:maximum(difference)

    hist!(
        ax,
        difference;
        color = (:purple, 0.5),
        bins = bins,
    )

    vlines!(
        ax,
        mean(difference);
        color = :black,
        linestyle = :dash,
        linewidth = 4
    )

    return fig
end

function singlescenario_test_positivity_plot(
    test_positivity_struct_vec; agg = :seven_day
)
    posoddsmatrix = reduce(
        hcat,
        map(array -> array[:, 1], getfield.(test_positivity_struct_vec, agg)),
    )
    avgpositivity = vec(mapslices(NaNMath.mean, posoddsmatrix; dims = 2))

    fig = Figure()
    ax = Axis(
        fig[1, 1]; xlabel = "Time steps by $(agg)", ylabel = "Test Positivity"
    )
    lines!(ax, 1:length(avgpositivity), avgpositivity)

    return fig
end

function test_positivity_distribution_plot(
    test_positivity_struct_vec; agg = :seven_day, kwargs...
)
    df = @chain test_positivity_struct_vec begin
        getfield.(agg)
        reduce(vcat, _)
        DataFrame(Tables.table(_), [:positivity, :outbreak])
    end

    df[!, :outbreak] = string.(df[:, :outbreak])

    return draw(
        data(df) *
        mapping(
            :positivity => "Test Positivity aggregated by $(agg)"; kwargs...
        ) *
        histogram(; bins = 0.0:0.05:1.05), ;
        axis = (ylabel = "Count",),
    )
end

function save_compare_ensemble_OTchars_plot(
    char_struct_vec,
    columnfacetchar::Symbol,
    plottingchars;
    plotname,
    plotsrootdir = plotsdir("ensemble/testing-comparison"),
    clinic_tested_dir,
    plotformat = "png",
    resolution = (2200, 1200),
    columnfacetchar_label = "Alert Threshold",
    binwidth = 1.0,
    xlabel = "Alert Characteristic Value",
    ylabel = "Density",
    legend = true,
    legendlabel = "Outbreak Chacteristic",
    meanlines = false,
    meanlabels = false,
    normalization = :none,
    kwargs...,
)
    plot = compare_ensemble_OTchars_plots(
        char_struct_vec,
        columnfacetchar,
        plottingchars;
        columnfacetchar_label = columnfacetchar_label,
        binwidth = binwidth,
        xlabel = xlabel,
        ylabel = ylabel,
        legend = legend,
        legendlabel = legendlabel,
        meanlines = meanlines,
        meanlabels = meanlabels,
        normalization = normalization,
        kwargs...,
    )

    plotpath = joinpath(
        plotsrootdir, clinic_tested_dir
    )
    mkpath(plotpath)

    save(
        joinpath(plotpath, "$plotname.$plotformat"),
        plot;
        resolution = resolution,
    )

    return nothing
end

function compare_ensemble_OTchars_plots(
    char_struct_vec,
    columnfacetchar::Symbol,
    plottingchars;
    columnfacetchar_label = "Alert Threshold",
    binwidth = 1.0,
    xlabel = "Alert Characteristic Value",
    ylabel = "Density",
    legend = true,
    legendlabel = "Outbreak Chacteristic",
    meanlines = false,
    meanlabels = false,
    normalization = :none,
    kwargs...,
)
    xs, ys = calculate_comparison_plot_facet_dims(
        char_struct_vec, columnfacetchar
    )
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    @pack! kwargs_dict = columnfacetchar_label,
    binwidth,
    xlabel,
    ylabel,
    legendlabel

    fig = Figure()

    construct_OTchars_facets!(
        fig,
        char_struct_vec,
        plottingchars,
        xs,
        ys,
        columnfacetchar,
        kwargs_dict;
        meanlines = meanlines,
        meanlabels = meanlabels,
        normalization = normalization,
    )

    if legend
        Legend(
            fig[:, end + 1],
            [
                PolyElement(; color = col) for
                col in map(chartuple -> chartuple.color, plottingchars)
            ],
            map(chartuple -> chartuple.label, plottingchars);
            label = legendlabel,
        )
    end
    return fig
end

function calculate_comparison_plot_facet_dims(
    char_struct_vec, facetchar
)
    xlength = length(
        Set(
            getfield.(
                getfield.(char_struct_vec, :outbreak_detect_spec),
                facetchar
            ),
        ),
    )
    ylength = length(
        Set(getfield.(getfield.(char_struct_vec, :ind_test_spec), :specificity))
    )

    xs = repeat(1:ylength, xlength)
    ys = repeat(1:xlength; inner = ylength)

    return xs, ys
end

function construct_OTchars_facets!(
    fig,
    char_struct_vec,
    plottingchars,
    xs,
    ys,
    columnfacetchar,
    kwargs_dict;
    meanlines = false,
    meanlabels = false,
    normalization = :pdf,
)
    for (OT_char_tuple, x, y) in zip(char_struct_vec, xs, ys)
        charvecs = map(
            chartuple -> reduce(
                vcat, getproperty(OT_char_tuple.OT_chars, chartuple.char)
            ),
            plottingchars,
        )

        @unpack binwidth, xlabel, ylabel, columnfacetchar_label = kwargs_dict

        if !haskey(kwargs_dict, :bins)
            bins = calculate_bins(charvecs, binwidth)
        else
            bins = kwargs_dict[:bins]
        end

        gl = fig[x, y] = GridLayout()
        ax = Axis(
            gl[2, 1];
            xlabel = xlabel,
            ylabel = ylabel,
        )

        for charnumber in eachindex(plottingchars)
            if isempty(charvecs[charnumber])
                break
            end
            hist!(
                ax,
                charvecs[charnumber];
                bins = bins,
                color = plottingchars[charnumber].color,
                normalization = normalization,
            )

            if meanlines || meanlabels
                charmean = mean(charvecs[charnumber])
            end
            if meanlines
                vlines!(
                    ax,
                    charmean;
                    color = :black,
                    linestyle = :dash,
                    linewidth = 4,
                )
            end
            if meanlabels
                hjust = 0
                vjust = 0
                if haskey(plottingchars[charnumber], :hjust)
                    hjust = plottingchars[charnumber].hjust
                end
                if haskey(plottingchars[charnumber], :vjust)
                    vjust = plottingchars[charnumber].vjust
                end
                text!(
                    Point(charmean + hjust, 0 + vjust);
                    text = "Mean ($(plottingchars[charnumber].label)):\n$(round(charmean, digits = 2))",
                )
            end
        end

        Label(
            gl[1, :],
            L"\text{\textbf{Individual Test} - Sensitivity: %$(OT_char_tuple.ind_test_spec.sensitivity), Specificity: %$(OT_char_tuple.ind_test_spec.specificity), %$(columnfacetchar_label): %$(getfield(OT_char_tuple.outbreak_detect_spec, columnfacetchar)), Perc Clinic Tested: %$(OT_char_tuple.outbreak_detect_spec.percent_clinic_tested)}";
            word_wrap = true,
        )
        colsize!(gl, 1, Relative(1))
    end
end

function calculate_bins(charvec, binwidth)
    filteredcharvec = filter(
        !isempty, charvec
    )
    minbinvec = minimum.(filteredcharvec)
    maxbinvec = maximum.(filteredcharvec)
    minbin = minimum(minbinvec)
    maxbin = maximum(maxbinvec)
    minbin -= 3 * binwidth / 2
    maxbin += 3 * binwidth / 3
    if minbin == maxbin
        minbin -= binwidth
        maxbin += binwidth
    end
    return minbin:binwidth:maxbin
end

function compare_optimal_thresholds_chars_plot(
    optimal_thresholds_vec,
    plottingchars;
    kwargs...
)
    unique_percent_clinic_tested = unique(
        optimal_thresholds_vec.percent_clinic_tested
    )

    clinical_case_test_spec = IndividualTestSpecification(1.0, 0.0)

    for percent_clinic_tested in unique_percent_clinic_tested
        optimal_thresholds_chars = filter(
            optimal_thresholds ->
                optimal_thresholds.percent_clinic_tested ==
                percent_clinic_tested ||
                    (
                        optimal_thresholds.percent_clinic_tested ==
                        1.0 &&
                        optimal_thresholds.individual_test_specification ==
                        clinical_case_test_spec
                    ),
            optimal_thresholds_vec,
        )

        plot = create_optimal_thresholds_chars_plot(
            optimal_thresholds_chars,
            plottingchars;
            kwargs...
        )

        plotpath = plotsdir(
            "ensemble/testing-comparison/clinic-tested_$percent_clinic_tested"
        )
        mkpath(plotpath)

        save(
            joinpath(
                plotpath,
                "compare-outreak_clinic-tested-$(percent_clinic_tested)_best-thresholds.png",
            ),
            plot;
            resolution = (2200, 1600),
        )

        @info "Created optimal thresholds plot for % clinic tested $(percent_clinic_tested)"
    end

    return nothing
end

function create_optimal_thresholds_chars_plot(
    optimal_thresholds_chars,
    plottingchars;
    kwargs...
)
    number_tests = length(optimal_thresholds_chars)

    sort!(
        optimal_thresholds_chars;
        by = threshold ->
            threshold.individual_test_specification.specificity,
    )

    fig = Figure()

    thresholdschars_structarr =
        optimal_thresholds_chars.outbreak_threshold_chars
    for (x, chartuple) in pairs(plottingchars)
        bins_vec = Vector{StepRangeLen}(undef, number_tests)
        thresholdschars_vec =
            getproperty.(thresholdschars_structarr, chartuple.char)

        charvecs = reduce(vcat, thresholdschars_vec)

        if !haskey(chartuple, :bins)
            if !haskey(chartuple, :binwidth)
                @error "The metric $(chartuple.char) wasn't provided with bins or a binwidth"
                break
            end
            bins = calculate_bins(charvecs, chartuple.binwidth)
        else
            bins .= chartuple.bins
        end

        if !haskey(chartuple, :label)
            label = :none
        else
            label = chartuple.label
        end

        for (y, optimal_thresholds) in pairs(optimal_thresholds_chars)
            gl = fig[y, x] = GridLayout()
            ax = Axis(gl[2, 1]; xlabel = label)

            hist!(
                ax,
                reduce(vcat, thresholdschars_vec[y]);
                bins = bins,
                color = chartuple.color,
            )
            Label(
                gl[1, :],
                L"\text{\textbf{Individual Test} - Sensitivity: %$(optimal_thresholds.individual_test_specification.sensitivity), Specificity: %$(optimal_thresholds.individual_test_specification.specificity), Alert Threshold: %$(optimal_thresholds.alert_threshold), Perc Clinic Tested: %$(optimal_thresholds.percent_clinic_tested)}";
                word_wrap = true,
            )
            colsize!(gl, 1, Relative(1))

            if y < number_tests
                hidexdecorations!(ax; ticklabels = false, ticks = false)
            end
        end
    end

    return fig
end
# end
#
