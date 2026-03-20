#%%
using DrWatson: DrWatson
using OutbreakDetectionCore: OutbreakDetectionCore
using OutbreakDetection
using Try: Try
using DataFrames
using CSV: CSV

# Make sure these values are present in the optimization script
alert_method = OutbreakDetectionCore.AlertMethod(OutbreakDetectionCore.MovingAverage(7))
accuracy_metric = OutbreakDetectionCore.AccuracyMetric(OutbreakDetectionCore.BalancedAccuracy())
threshold_bounds = (; lower = 0.0, upper = 50.0)
alert_filtering_strategy = OutbreakDetectionCore.AlertFilteringStrategy(OutbreakDetectionCore.AllAlerts())
alert_outbreak_matching_strategy = OutbreakDetectionCore.AlertOutbreakMatchingStrategy(OutbreakDetectionCore.SingleOutbreakPerAlert())

mkpath(supplement_tables())

#%%
optimized_filedir = OutbreakDetectionCore.outdir("ensemble", "threshold-optimization")
optimized_threshold_results = Try.unwrap(
    OutbreakDetectionCore.load_previous_optimization_results_structvector(
        optimized_filedir,
        "threshold-optimization.jld2",
        joinpath(optimized_filedir, "checkpoints");
        default_action = :continue
    )
);

#%%
filtered_results = filter(
    r -> r.alert_method == alert_method &&
        r.accuracy_metric == accuracy_metric &&
        r.threshold_bounds == threshold_bounds &&
        r.alert_filtering_strategy == alert_filtering_strategy &&
        r.alert_outbreak_matching_strategy == alert_outbreak_matching_strategy,
    optimized_threshold_results
);

#%%
# Create the wide table directly from the vector of optimization results.
wide_thresholds_df = create_wide_optimal_thresholds_df(
    filtered_results,
    :optimal_threshold;
    simplify = true,
)

#%%
unique_noise_levels = sort(unique(optimized_threshold_results.noise_level))

function mean_delay_difference_vs_perfect_test_0_day(df::DataFrame)
    percent_columns = filter(
        column_name -> try
            parse(Float64, string(column_name))
            return true
        catch
            return false
        end,
        names(df),
    )

    baseline_row = only(
        eachrow(
            DataFrames.subset(
                df,
                :test_type => ByRow(==("Perfect Test")),
                :test_lag => ByRow(==(0)),
            ),
        ),
    )

    baseline_values = [parse(Float64, string(baseline_row[column])) for column in percent_columns]

    summary_df = DataFrames.subset(
        df,
        [:test_type, :test_lag] =>
            ByRow((test_type, test_lag) -> !(test_type == "Perfect Test" && test_lag == 0)),
    )

    summary_df[!, :mean_delay_difference_vs_perfect_test_0_day] = map(
        eachrow(summary_df),
    ) do row
        row_values = [parse(Float64, string(row[column])) for column in percent_columns]
        return sum(row_values .- baseline_values) / length(percent_columns)
    end

    return DataFrames.select(
        summary_df,
        :noise_type,
        :noise_level,
        :test_type,
        :test_lag,
        :mean_delay_difference_vs_perfect_test_0_day,
    )
end

function prepare_delay_difference_summary_format(df::DataFrame)
    formatted_df = copy(df)
    formatted_df.noise_type .= String.(formatted_df.noise_type)
    replace!(
        formatted_df.noise_type,
        "static" => "Static Noise",
        "dynamic" => "Dynamic Noise",
        "all_noise_structures" => "All noise structures",
    )

    DataFrames.rename!(
        formatted_df,
        :noise_type => "Noise Type",
        :test_type => "Test Type",
        :test_lag => "Test Lag",
        :mean_delay_difference_vs_perfect_test_0_day =>
            "Mean Delay Difference vs Perfect Test (0-day lag)",
    )

    return formatted_df
end

#%%
for noise_level in unique_noise_levels

    local df = subset_for_noise_level_with_perfect_tests(wide_thresholds_df, noise_level)
    cleaned_df = prepare_wide_optimal_thresholds_df_format(df)
    Base.display(cleaned_df)

    if noise_level == 7.0
        CSV.write(
            supplement_tables("optimal-thresholds.csv"),
            DataFrames.select!(cleaned_df, DataFrames.Not(:noise_level))
        )
    end

end

#%%
wide_accuracy_df = create_wide_optimal_thresholds_df(
    filtered_results,
    :accuracies;
    simplify = true,
)

#%%
for noise_level in unique_noise_levels

    local df = subset_for_noise_level_with_perfect_tests(wide_accuracy_df, noise_level)
    cleaned_df = prepare_wide_optimal_thresholds_df_format(df)
    Base.display(cleaned_df)
end

#%%
wide_delay_df = create_wide_optimal_thresholds_df(
    filtered_results,
    :detection_delays;
    simplify = true,
)

#%%
for noise_level in unique_noise_levels
    local df = subset_for_noise_level_with_perfect_tests(wide_delay_df, noise_level)
    Base.display(prepare_wide_optimal_thresholds_df_format(df))

    cleaned_df = prepare_delay_difference_summary_format(
        mean_delay_difference_vs_perfect_test_0_day(df),
    )
    Base.display(cleaned_df)
end

#%%
filtered_multi_results = filter(
    r -> r.alert_method == alert_method &&
        r.accuracy_metric == accuracy_metric &&
        r.threshold_bounds == threshold_bounds &&
        r.alert_filtering_strategy == alert_filtering_strategy &&
        r.alert_outbreak_matching_strategy == OutbreakDetectionCore.AlertOutbreakMatchingStrategy(OutbreakDetectionCore.MultipleOutbreaksPerAlert()),
    optimized_threshold_results
);

#%%
# Create the wide table directly from the vector of optimization results.
multi_wide_thresholds_df = create_wide_optimal_thresholds_df(
    filtered_multi_results,
    :optimal_threshold;
    simplify = true,
)

#%%
for noise_level in unique_noise_levels
    local df = subset_for_noise_level_with_perfect_tests(multi_wide_thresholds_df, noise_level)
    cleaned_df = prepare_wide_optimal_thresholds_df_format(df)
    Base.display(cleaned_df)
end

#%%
multi_wide_accuracy_df = create_wide_optimal_thresholds_df(
    filtered_multi_results,
    :accuracies;
    simplify = true,
)

#%%
for noise_level in unique_noise_levels
    local df = subset_for_noise_level_with_perfect_tests(multi_wide_accuracy_df, noise_level)
    cleaned_df = prepare_wide_optimal_thresholds_df_format(df)
    Base.display(cleaned_df)
end
