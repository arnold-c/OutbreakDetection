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
threshold_bounds = (; lower = 0.0, upper = 20.0)
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

wide_accuracy_df = create_wide_optimal_thresholds_df(
    filtered_results,
    :accuracies;
    simplify = true,
)


#%%
unique_noise_levels = unique(optimized_threshold_results.noise_level)

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
