using DataFrames

function line_accuracy_plot(
    optimal_thresholds_vec;
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

    println(long_df)

    unique_tests = select(
        unique(long_df, [:sensitivity, :test_lag]), [:sensitivity, :test_lag]
    )

    println(unique_tests)

    fig = Figure()

    _line_accuracy_facet!(
        fig,
        long_df;
        colors = colors,
    )

    Legend(
        fig[0, :],
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
    rowsize!(fig.layout, 0, Relative(0.1))
    colsize!(fig.layout, 1, Relative(1))

    return fig
end

function _line_accuracy_facet!(
    fig,
    long_df;
    colors = Makie.wong_colors(),
)
    # gl = fig[y + 2, x + 1] = GridLayout()
    gl = fig[1, 1] = GridLayout()
    ax = Axis(
        gl[1, 1];
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

    return nothing
end
