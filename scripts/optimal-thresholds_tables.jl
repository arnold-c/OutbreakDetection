#%%
using DrWatson: DrWatson
using OutbreakDetectionCore: OutbreakDetectionCore
using OutbreakDetection
using Try: Try
using DataFrames

# Make sure these values are present in the optimization script
alert_method = OutbreakDetectionCore.AlertMethod(OutbreakDetectionCore.MovingAverage(7))
accuracy_metric = OutbreakDetectionCore.AccuracyMetric(OutbreakDetectionCore.BalancedAccuracy())
threshold_bounds = (; lower = 0.0, upper = 20.0)
alert_filtering_strategy = OutbreakDetectionCore.AlertFilteringStrategy(OutbreakDetectionCore.AllAlerts())
plotdirpath = DrWatson.plotsdir()

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
        r.alert_filtering_strategy == alert_filtering_strategy,
    optimized_threshold_results
)
# Reshape results into matrix structure
result_matrix, unique_noise_types = OutbreakDetection.reshape_optimization_results_to_matrix(filtered_results);

unique_test_specifications = OutbreakDetection.get_unique_test_specifications_in_sorted_order(
    result_matrix[1, 1].test_specification
)

unique_noise_levels = unique(optimized_threshold_results.noise_level)

#%%
function calculate_optimal_threshold_summaries(
        char_vecs; percentiles = [0.25, 0.5, 0.75], nboots = 10000, ci = 0.95
    )
    all_chars = reduce(vcat, char_vecs)

    if isempty(all_chars)
        return missing, repeat([missing], length(percentiles))
    end

    if !isnothing(percentiles)
        char_percentiles = map(
            percentile -> StatsBase.quantile(all_chars, percentile), percentiles
        )
    else
        char_percentiles = [missing]
    end

    if !isnothing(nboots) && !isnothing(ci)
        char_mean, char_lower, char_upper = Bootstrap.confint(
            Bootstrap.bootstrap(
                StatsBase.mean, all_chars, Bootstrap.BasicSampling(nboots)
            ),
            Bootstrap.BCaConfInt(ci),
        )[1]
    else
        char_mean = StatsBase.mean(all_chars)
        char_lower = missing
        char_upper = missing
    end

    return char_lower, char_mean, char_upper, char_percentiles
end

#%%
function create_optimal_threshold_summary_df(
        optimal_thresholds_vec,
        characteristics::T;
        percentiles = [0.25, 0.5, 0.75],
        nboots = 10000,
        ci = 0.95,
    ) where {T <: AbstractVector{<:Symbol}}
    if !isnothing(percentiles)
        percentile_labels = map(percentiles) do perc
            perc = perc * 100
            try
                return "$(Int64(perc))th"
            catch p
                return replace("$(perc)th", "." => "_")
            end
        end
    else
        percentile_labels = ["100th"]
    end

    chars_df = mapreduce(hcat, characteristics) do characteristic
        char_percentile_labels = map(
            perc -> "$(characteristic)_$perc", percentile_labels
        )
        char_mean_label = "$(characteristic)_mean"

        ci_label = "CI_$(Int64(ci * 100))"
        char_mean_lower_label = "$(characteristic)_lower_$(ci_label)"
        char_mean_upper_label = "$(characteristic)_upper_$(ci_label)"

        DataFrames.DataFrame(
            mapreduce(vcat, optimal_thresholds_vec) do opt
                char_lower, char_mean, char_upper, chars_percentiles = calculate_optimal_threshold_summaries(
                    getfield(opt, characteristic);
                    percentiles = percentiles,
                    nboots = nboots,
                    ci = ci,
                )

                return char_lower,
                    char_mean, char_upper,
                    chars_percentiles...
            end,
            [
                char_mean_lower_label,
                char_mean_label,
                char_mean_upper_label,
                char_percentile_labels...,
            ],
        )
    end

    if isnothing(nboots) || isnothing(ci)
        DataFrames.select!(chars_df, DataFrames.Not(r"CI"))
    end

    if isnothing(percentiles)
        DataFrames.select!(chars_df, DataFrames.Not(r".*_[0-9]+th"))
    end

    if isnothing(percentiles) && (isnothing(nboots) || isnothing(ci))
        DataFrames.rename!(
            x -> replace(x, "_mean" => ""),
            chars_df;
            cols = DataFrames.Cols(r".*_mean"),
        )
    end

    core_df = DataFrames.DataFrame(
        "noise_type" =>
            getfield.(optimal_thresholds_vec, :noise_type_description),
        "noise_level" =>
            getfield.(optimal_thresholds_vec, :noise_level),
        "percent_tested" =>
            getfield.(optimal_thresholds_vec, :percent_tested),
        "sensitivity" =>
            getfield.(
            optimal_thresholds_vec.test_specification,
            :sensitivity,
        ),
        "specificity" =>
            getfield.(
            optimal_thresholds_vec.test_specification,
            :specificity,
        ),
        "test_lag" =>
            getfield.(
            optimal_thresholds_vec.test_specification,
            :test_result_lag,
        ),
        "alert_threshold" =>
            getfield.(
            optimal_thresholds_vec,
            :optimal_threshold,
        ),
        "accuracy" => getfield.(
            optimal_thresholds_vec,
            :accuracies,
        ),
    )

    return hcat(core_df, chars_df)
end

#%%
thresholds_summary_df = create_optimal_threshold_summary_df(
    optimized_threshold_results,
    [:unavoidable_cases, :detection_delays, :alert_durations];
    percentiles = nothing,
    nboots = nothing,
)

#%%
wide_thresholds_df = OutbreakDetectionCore.create_wide_optimal_thresholds_df(
    thresholds_summary_df,
    :alert_threshold
)

#%%
DataFrames.subset(
    wide_thresholds_df,
    [:sensitivity, :specificity] => (sens, spec) -> sens .== 1.0 .&& spec .== 1.0
) |>
    df -> for (i, col) in pairs(names(df)[6:end])
    if length(unique(df[!, col])) != 1
        @show col
        @show df[!, Cols(1:5, col)]
        @show unique(df[!, col])
    end
end

#%%
DataFrames.subset(
    wide_thresholds_df,
    [:noise_level, :sensitivity, :specificity, :test_lag] => (noise, sens, spec, lag) -> noise .== 1.0 .&& sens .== 1.0 .&& spec .== 1.0 .&& lag .== 0
)

# This is what is needs to match
#| All noise | Perfect | 0 | 1.172 | 2.734 | 3.516 | 5.859 | 6.641 | 7.422 | 8.984 | 10.547 | 11.328 | 12.109 |
