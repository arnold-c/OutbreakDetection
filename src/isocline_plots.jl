using LinearAlgebra

function isocline_accuracy_plot(
    optimal_thresholds_vec,
    ::T,
    isocline = 0.85;
    text_color_threshold = 0.9,
    digits = 2,
) where {T<:Val{Heatmap}}
    long_df = create_optimal_thresholds_df(optimal_thresholds_vec)

    unique_testing_rates = sort(unique(long_df.percent_clinic_tested))

    filtered_df = DataFrames.subset(long_df, :accuracy => x -> x .>= isocline)

    wide_df = create_wide_optimal_thresholds_df(filtered_df, :accuracy)
    val_mat =
        round.(
            Transpose(
                Matrix(
                    select(
                        wide_df, Not([:sensitivity, :specificity, :test_lag])
                    ),
                ),
            );
            digits = digits,
        )
    replace!(val_mat, missing => NaN)
    println(wide_df)

    fig = Figure()
    ax = Axis(
        fig[1, 1]; xlabel = "Testing Rate", ylabel = "Test Characteristic"
    )

    for test in axes(val_mat, 2)
        min_rate = 1
        min_rate_val = val_mat[min_rate, test]
        for testing_rate in axes(val_mat, 1)
            val = val_mat[testing_rate, test]
            if isnan(min_rate_val)
                min_rate = testing_rate
                min_rate_val = val
            end
            if testing_rate > min_rate
                val_mat[testing_rate, test] = NaN
            end
            filtered_val = val_mat[testing_rate, test]
            if isnan(filtered_val)
                continue
            end
        end
    end

    heatmap!(
        ax,
        val_mat;
        colormap = :Greens,
    )

    for test in axes(val_mat, 2)
        for testing_rate in axes(val_mat, 1)
            val = val_mat[testing_rate, test]
            if isnan(val)
                continue
            end
            txtcolor = val >= text_color_threshold ? :white : :black
            text!(
                ax,
                # unique_testing_rates[testing_rate],
                # collect(1:6)[test],
                "$(val)";
                position = (testing_rate, test),
                color = txtcolor,
                align = (:center, :center),
            )
        end
    end

    return fig
end

function isocline_accuracy_plot(
    optimal_thresholds_vec,
    ::T,
    isocline = 0.85;
    digits = 2,
) where {T<:Val{Scatter}}
    long_df = create_optimal_thresholds_df(optimal_thresholds_vec)

    filtered_df = DataFrames.subset(long_df, :accuracy => x -> x .>= isocline)

    wide_df = create_wide_optimal_thresholds_df(filtered_df, :accuracy)
    val_mat =
        round.(
            Transpose(
                Matrix(
                    select(
                        wide_df, Not([:sensitivity, :specificity, :test_lag])
                    ),
                ),
            );
            digits = digits,
        )
    replace!(val_mat, missing => NaN)
    println(wide_df)

    fig = Figure()
    ax = Axis(
        fig[1, 1]; xlabel = "Testing Rate", ylabel = "Test Characteristic"
    )

    unique_percent_clinic_tested = sort(
        unique(filtered_df.percent_clinic_tested)
    )
    testing_rate_vec = zeros(Float64, length(axes(val_mat, 2)))
    for test in axes(val_mat, 2)
        for testing_rate in axes(val_mat, 1)
            if isnan(val_mat[testing_rate, test])
                continue
            else
                testing_rate_vec[test] = unique_percent_clinic_tested[testing_rate]
                break
            end
        end
    end

    lines!(
        ax,
        testing_rate_vec,
        axes(val_mat, 2);
        color = :black,
        linewidth = 3,
    )

    return fig
end
