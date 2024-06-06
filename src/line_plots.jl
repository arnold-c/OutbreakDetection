using DataFrames

function line_accuracy_plot(
    optimal_thresholds_vec
)
    long_df = create_optimal_threshold_summary_df(
        optimal_thresholds_vec,
        :accuracy;
        percentiles = [0.1, 0.5, 0.9],
    )

    select!(long_df, Cols(
        :percent_clinic_tested,
        :sensitivity,
        :test_lag,
        x -> endswith(x, "th"),
        )
    )

    println(long_df)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "Testing Rate", ylabel = "Accuracy")

    colors = Makie.wong_colors()

    for (i, test) in pairs(unique(long_df.sensitivity))
        subsetted_df = DataFrames.subset(long_df, :sensitivity => x -> x .== test)

        lines!(
            ax,
            subsetted_df.percent_clinic_tested,
            subsetted_df.accuracy_50th;
            color = colors[i],
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

    return fig
end
