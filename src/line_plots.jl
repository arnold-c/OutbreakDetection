using DataFrames

function line_accuracy_plot(
    noise_spec_vec,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params,
)
    fig = Figure()

    noise_descriptions = get_noise_description.(noise_spec_vec)
    unique_noise_descriptions = unique(noise_descriptions)

    for (j, noise_description) in pairs(unique_noise_descriptions)
        shape_noise_specification = filter(
            noise_spec ->
                noise_description == get_noise_description(noise_spec),
            noise_spec_vec,
        )

        @show j
        @show noise_description
        # @show shape_noise_specification

        for (i, noise_spec) in pairs(shape_noise_specification)
            @show i, noise_spec
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
                optimal_thresholds_vec,
                fig,
                j,
                i,
            )
        end
    end

    return fig
end

function _line_accuracy_plot!(
    optimal_thresholds_vec,
    fig,
    j,
    i;
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

    unique_tests = select(
        unique(long_df, [:sensitivity, :test_lag]), [:sensitivity, :test_lag]
    )

    noise_spec = optimal_thresholds_vec[1].noise_specification

    gl = fig[j, i] = GridLayout(3, 2)

    _line_accuracy_facet!(
        gl,
        noise_spec,
        long_df;
        colors = colors,
    )

    Legend(
        gl[1, :],
        [
            PolyElement(; color = colors[i]) for
            i in eachindex(unique_tests.sensitivity)
        ],
        map(
            test -> "$(test.sensitivity), $(test.test_lag)",
            eachrow(unique_tests),
        ),
        "Test Type";
        orientation = :horizontal,
    )
    rowsize!(gl, 1, Relative(0.1))

    return nothing
end

function _line_accuracy_facet!(
    gl,
    noise_spec,
    long_df;
    colors = Makie.wong_colors(),
)
    ax = Axis(
        gl[3, 2];
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

    Box(gl[2, 2]; color = :lightgray, strokevisible = false)
    Label(
        gl[2, 2],
        "Noise type: $(noise_spec.noise_type)";
        fontsize = 20,
        padding = (0, 0, 5, 5),
        valign = :bottom,
    )

    Box(
        gl[3, 1];
        color = :lightgray,
        strokevisible = false,
    )
    Label(
        gl[3, 1],
        "Noise scaling: $(noise_spec.noise_mean_scaling)";
        fontsize = 20,
        rotation = pi / 2,
        padding = (5, 5, 0, 0),
        valign = :center,
    )

    rowsize!(gl, 3, Relative(0.85))
    colsize!(gl, 2, Relative(0.97))
    return nothing
end
