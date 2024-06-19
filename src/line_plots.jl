using DataFrames
using DrWatson: DrWatson
using StatsBase: StatsBase

lineplot_colors = [
    "#56B4E9"
    "#E69F00"
    repeat(["#000000"], 2)...
]

function line_accuracy_plot(
    noise_spec_vec,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    plotdirpath = DrWatson.plotsdir(),
    plotname = "line_accuracy_plot",
    plotformat = "png",
    size = (2200, 1200),
    colors = lineplot_colors,
    xlabel = "Proportion Tested",
    ylabel = "Accuracy",
    labelsize = 24,
    show_x_facet_label = true,
    show_y_facet_label = true,
    force = false,
    kwargs...,
)
    mkpath(plotdirpath)
    plotpath = joinpath(plotdirpath, "$plotname.$plotformat")

    if !isfile(plotpath) || force
        fig = Figure()

        noise_descriptions = get_noise_description.(noise_spec_vec)
        unique_noise_descriptions = unique(noise_descriptions)
        num_noise_descriptions = length(unique_noise_descriptions)
        unique_test_specifications = unique(optimal_threshold_test_spec_vec)

        for (i, noise_description) in pairs(unique_noise_descriptions)
            shape_noise_specification = filter(
                noise_spec ->
                    noise_description == get_noise_description(noise_spec),
                noise_spec_vec,
            )

            if contains(noise_description, "dynamical")
                shape_noise_specification = reverse(shape_noise_specification)
            end

            for (j, noise_spec) in pairs(shape_noise_specification)
                optimal_threshold_comparison_params = (
                    noise_specification = noise_spec,
                    optimal_threshold_core_params...,
                )

                optimal_thresholds_vec = calculate_OptimalThresholdCharacteristics(
                    ensemble_percent_clinic_tested_vec,
                    optimal_threshold_test_spec_vec,
                    optimal_threshold_comparison_params,
                )

                _line_accuracy_plot!(
                    fig,
                    noise_spec,
                    unique_test_specifications,
                    optimal_thresholds_vec,
                    i,
                    j;
                    num_noise_descriptions = num_noise_descriptions,
                    colors = colors,
                    xlabel = xlabel,
                    ylabel = ylabel,
                    show_x_facet_label = show_x_facet_label,
                    kwargs...,
                )
            end
        end

        if show_y_facet_label
            map(enumerate(unique_noise_descriptions)) do (i, noise_description)
                Box(fig[i, 0]; color = :lightgray, strokevisible = false)
                Label(
                    fig[i, 0],
                    titlecase(noise_description);
                    fontsize = 16,
                    rotation = pi / 2,
                    padding = (0, 0, 0, 0),
                    valign = :center,
                    tellheight = false,
                )
            end
            colsize!(fig.layout, 0, Relative(0.03))
        end

        Legend(
            fig[0, :],
            map(enumerate(unique_test_specifications)) do (i, test_spec)
                linestyle = test_spec.test_result_lag == 0 ? :solid : :dot
                return [
                    PolyElement(; color = (colors[i], 0.3)),
                    LineElement(; color = colors[i], linestyle = linestyle),
                ]
            end,
            get_test_description.(unique_test_specifications);
            labelsize = labelsize,
            orientation = :horizontal,
        )
        rowsize!(fig.layout, 0, Relative(0.03))

        Makie.save(plotpath, fig; size = size)
        return fig
    end

    return nothing
end

function _line_accuracy_plot!(
    fig,
    noise_spec,
    unique_test_specifications,
    optimal_thresholds_vec,
    i,
    j;
    colors = lineplot_colors,
    xlabel = "Proportion Tested",
    ylabel = "Accuracy",
    show_x_facet_label = true,
    num_noise_descriptions = 1,
    kwargs...,
)
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    long_df = create_optimal_threshold_summary_df(
        optimal_thresholds_vec,
        :accuracy;
        percentiles = [0.1, 0.9],
    )

    if show_x_facet_label
        x_facet_label = "Mean daily noise: $(round(
            StatsBase.mean(optimal_thresholds_vec[1].outbreak_threshold_chars.mean_poisson_noise) /
            StatsBase.mean(optimal_thresholds_vec[1].outbreak_threshold_chars.poisson_noise_prop);
            digits = 2,
        ))"

        kwargs_dict[:x_facet_label] = x_facet_label
    end

    select!(
        long_df,
        Cols(
            :percent_clinic_tested,
            :sensitivity,
            :specificity,
            :test_lag,
            :accuracy_mean,
            x -> endswith(x, "th"),
        ),
    )

    gl = fig[i, j] = GridLayout()

    if i != num_noise_descriptions
        xlabel = ""
    end

    if j != 1
        ylabel = ""
    end

    _line_accuracy_facet!(
        gl,
        noise_spec,
        unique_test_specifications,
        long_df;
        colors = colors,
        xlabel = xlabel,
        ylabel = ylabel,
        kwargs_dict...,
    )

    return nothing
end

function _line_accuracy_facet!(
    gl,
    noise_spec,
    unique_test_specifications,
    long_df;
    colors = lineplot_colors,
    xlabel = "Proportion Tested",
    ylabel = "Accuracy",
    kwargs...,
)
    kwargs_dict = Dict(kwargs)

    ypos = haskey(kwargs_dict, :x_facet_label) ? 2 : 1
    xpos = 1

    ax = Axis(
        gl[ypos, xpos];
        xlabel = xlabel,
        ylabel = ylabel,
    )

    for (i, test) in pairs(unique_test_specifications)
        subsetted_df = DataFrames.subset(
            long_df,
            :sensitivity =>
                x -> x .== test.sensitivity,
            :specificity =>
                x -> x .== test.specificity,
            :test_lag => x -> x .== test.test_result_lag,
        )

        linestyle = test.test_result_lag == 0 ? :solid : :dash

        string(
            subsetted_df.sensitivity[1], ", ", subsetted_df.test_lag[1]
        )
        lines!(
            ax,
            subsetted_df.percent_clinic_tested,
            subsetted_df.accuracy_mean;
            color = colors[i],
            linestyle = linestyle,
        )

        band!(
            ax,
            subsetted_df.percent_clinic_tested,
            subsetted_df.accuracy_10th,
            subsetted_df.accuracy_90th;
            color = colors[i],
            alpha = 0.3,
        )

        ylims!(ax, (0.6, 1))
    end

    if haskey(kwargs_dict, :x_facet_label)
        Box(gl[1, xpos]; color = :lightgray, strokevisible = false)
        Label(
            gl[1, xpos],
            kwargs_dict[:x_facet_label];
            fontsize = 16,
            padding = (0, 0, 0, 0),
            valign = :bottom,
            tellwidth = false,
        )
        rowsize!(gl, 2, Relative(0.9))
    end

    return nothing
end
