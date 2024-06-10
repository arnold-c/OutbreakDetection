using DataFrames
using DrWatson: DrWatson

function line_accuracy_plot(
    noise_spec_vec,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    plotdirpath = DrWatson.plotsdir(),
    plotname = "line_accuracy_plot",
    plotformat = "png",
    size = (2200, 1200),
    colors = Makie.wong_colors(),
    force = false,
)
    mkpath(plotdirpath)
    plotpath = joinpath(plotdirpath, "$plotname.$plotformat")

    if !isfile(plotpath) || force
        fig = Figure()

        noise_descriptions = get_noise_description.(noise_spec_vec)
        unique_noise_descriptions = unique(noise_descriptions)
        unique_test_descriptions =
            get_test_description.(unique(optimal_threshold_test_spec_vec))

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
                    optimal_thresholds_vec,
                    i,
                    j,
                )
            end
        end

        Legend(
            fig[0, :],
            map(
                i -> PolyElement(; color = colors[i]),
                eachindex(unique_test_descriptions),
            ),
            unique_test_descriptions,
            "Test Type";
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
    optimal_thresholds_vec,
    i,
    j;
    colors = Makie.wong_colors(),
)
    long_df = create_optimal_threshold_summary_df(
        optimal_thresholds_vec,
        :accuracy;
        percentiles = [0.1, 0.9],
    )

    select!(
        long_df,
        Cols(
            :percent_clinic_tested,
            :sensitivity,
            :test_lag,
            :accuracy_mean,
            x -> endswith(x, "th"),
        ),
    )

    gl = fig[i, j] = GridLayout()

    _line_accuracy_facet!(
        gl,
        noise_spec,
        long_df;
        colors = colors,
    )

    return nothing
end

function _line_accuracy_facet!(
    gl,
    noise_spec,
    long_df;
    colors = Makie.wong_colors(),
)
    ax = Axis(
        gl[2, 2];
        xlabel = "Testing Rate",
        ylabel = "Accuracy",
    )

    for (i, test) in pairs(unique(long_df.sensitivity))
        subsetted_df = DataFrames.subset(
            long_df, :sensitivity => x -> x .== test
        )
        test_label = string(
            subsetted_df.sensitivity[1], ", ", subsetted_df.test_lag[1]
        )
        lines!(
            ax,
            subsetted_df.percent_clinic_tested,
            subsetted_df.accuracy_mean;
            color = colors[i],
        )

        band!(
            ax,
            subsetted_df.percent_clinic_tested,
            subsetted_df.accuracy_10th,
            subsetted_df.accuracy_90th;
            color = colors[i],
            alpha = 0.3,
            label = test_label,
        )

        ylims!(ax, (0.6, 1))
    end

    Box(gl[1, 2]; color = :lightgray, strokevisible = false)
    Label(
        gl[1, 2],
        "Noise type: $(noise_spec.noise_type)";
        fontsize = 16,
        padding = (0, 0, 0, 0),
        valign = :bottom,
    )

    Box(
        gl[2, 1];
        color = :lightgray,
        strokevisible = false,
    )
    Label(
        gl[2, 1],
        get_noise_magnitude(noise_spec);
        fontsize = 16,
        rotation = pi / 2,
        padding = (0, 0, 0, 0),
        valign = :center,
    )

    rowsize!(gl, 2, Relative(0.9))
    colsize!(gl, 2, Relative(0.92))
    return nothing
end
