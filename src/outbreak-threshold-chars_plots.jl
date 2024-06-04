using GLMakie

function ensemble_OTChars_plot(
    OTChars,
    testspec,
    detectspec,
    plottingchars;
    plottitle = "",
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
    OT_char_tuple = (;
        OT_chars = OTChars,
        ind_test_spec = testspec,
        outbreak_detect_spec = detectspec,
    )
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    @pack! kwargs_dict = binwidth, xlabel, ylabel, legendlabel

    charvecs = map(
        chartuple -> reduce(
            vcat, getproperty(OT_char_tuple.OT_chars, chartuple.char)
        ),
        plottingchars,
    )

    fig = Figure()
    ax = Axis(fig[2, 1]; xlabel = xlabel, ylabel = ylabel)

    construct_single_OTchars_facet!(
        ax,
        charvecs,
        plottingchars,
        kwargs_dict;
        meanlines = meanlines,
        meanlabels = meanlabels,
        normalization = normalization,
    )

    if legend
        Legend(
            fig[2, 2],
            [
                PolyElement(; color = col) for
                col in map(chartuple -> chartuple.color, plottingchars)
            ],
            [chartuple.label for chartuple in plottingchars];
            label = legendlabel,
        )
    end

    Label(
        fig[1, :, Top()],
        plottitle,
    )

    rowsize!(fig.layout, 1, 5)

    return fig
end

function save_compare_ensemble_OTchars_plot(
    char_struct_vec,
    columnfacetchar::Symbol,
    plottingchars,
    percent_clinic_tested;
    plotname,
    plotdirpath,
    plotformat = "png",
    size = (2200, 1200),
    columnfacetchar_label = "Alert Threshold",
    binwidth = 1.0,
    xlabel = "Alert Characteristic Value",
    ylabel = "Density",
    legend = true,
    legendlabel = "Outbreak Chacteristic",
    meanlines = false,
    meanlabels = false,
    normalization = :none,
    force = false,
    kwargs...,
)
    mkpath(plotdirpath)

    plotpath = joinpath(plotdirpath, "$plotname.$plotformat")

    if !isfile(plotpath) || force
        plot = compare_ensemble_OTchars_plots(
            char_struct_vec,
            columnfacetchar,
            plottingchars,
            percent_clinic_tested;
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

        save(plotpath, plot; size = size)

        Makie.empty!(plot)
    end

    return nothing
end

function compare_ensemble_OTchars_plots(
    char_struct_vec,
    columnfacetchar::Symbol,
    plottingchars,
    percent_clinic_tested;
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

    Label(
        fig[1, :, Top()],
        "Perc Clinic Tested: $(percent_clinic_tested), Noise: $(getdirpath(char_struct_vec[end].noise_specification))",
    )

    unique_thresholds = unique(
        getfield.(
            getfield.(char_struct_vec, :outbreak_detect_spec), :alert_threshold
        ),
    )

    for (j, threshold) in pairs(unique_thresholds)
        Label(
            fig[2, j + 1, Top()],
            "Alert Threshold: $(threshold)",
        )
    end

    unique_test_specs = unique(getfield.(char_struct_vec, :ind_test_spec))

    for (i, test_spec) in pairs(unique_test_specs)
        Label(
            fig[i + 2, 1, Left()],
            "Sens: $(test_spec.sensitivity), Spec: $(test_spec.specificity),\nLag: $(test_spec.test_result_lag)";
            rotation = pi / 2,
        )
    end

    rowsize!(fig.layout, 1, 5)
    rowsize!(fig.layout, 2, 7)
    colsize!(fig.layout, 1, 7)

    return fig
end

function calculate_comparison_plot_facet_dims(
    char_struct_vec, facetchar
)
    xlength = length(
        Set(
            getfield.(
                getfield.(char_struct_vec, :outbreak_detect_spec),
                facetchar,
            ),
        ),
    )
    ylength = length(
        Set(getfield.(char_struct_vec, :ind_test_spec))
    )

    ys = repeat(1:ylength, xlength)
    xs = repeat(1:xlength; inner = ylength)

    return xs, ys
end

function construct_single_OTchars_facet!(
    ax,
    charvecs,
    plottingchars,
    kwargs_dict;
    meanlines = false,
    meanlabels = false,
    normalization = :pdf,
)
    @unpack binwidth = kwargs_dict

    if !haskey(kwargs_dict, :bins)
        bins = calculate_bins(charvecs, binwidth)
    else
        bins = kwargs_dict[:bins]
    end

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
end

function construct_OTchars_facets!(
    fig,
    char_struct_vec,
    plottingchars,
    xs,
    ys,
    kwargs_dict;
    meanlines = false,
    meanlabels = false,
    normalization = :pdf,
)
    number_tests = length(
        Set(getfield.(char_struct_vec, :ind_test_spec))
    )

    for (OT_char_tuple, x, y) in zip(char_struct_vec, xs, ys)
        skipped_plottingchar = 0

        charvecs = Vector{Vector{Union{Int64,Float64}}}(
            undef, length(plottingchars)
        )

        for (i, chartuple) in pairs(plottingchars)
            charvec = getproperty(OT_char_tuple.OT_chars, chartuple.char)

            if sum(isempty.(charvec)) == length(charvec)
                skipped_plottingchar += 1
                charvecs[i] = eltype(charvec)[]
            else
                charvecs[i] = reduce(vcat, charvec)
            end
        end

        @unpack xlabel, ylabel = kwargs_dict

        gl = fig[y + 2, x + 1] = GridLayout()
        ax = Axis(
            gl[1, 1];
            xlabel = xlabel,
            ylabel = ylabel,
        )

        if y < number_tests
            hidexdecorations!(ax; ticklabels = false, ticks = false)
        end

        hideydecorations!(ax)

        if skipped_plottingchar == length(plottingchars)
            continue
        end

        construct_single_OTchars_facet!(
            ax,
            charvecs,
            plottingchars,
            kwargs_dict;
            meanlines = meanlines,
            meanlabels = meanlabels,
            normalization = normalization,
        )
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
        linewidth = 4,
    )

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

