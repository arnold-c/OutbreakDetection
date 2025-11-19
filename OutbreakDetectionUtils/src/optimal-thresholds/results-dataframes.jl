export create_optimal_thresholds_df, create_wide_optimal_thresholds_df,
    create_optimal_threshold_summary_df, create_wide_optimal_threshold_summary_df,
    create_all_wide_optimal_threshold_summary_dfs

"""
    create_optimal_thresholds_df(optimal_thresholds_vec)

Create a long-format DataFrame from optimal threshold characteristics.
"""
function create_optimal_thresholds_df(optimal_thresholds_vec)
    return Chain.@chain begin
        map(optimal_thresholds_vec) do opt
            percent_clinic_tested = opt.percent_clinic_tested

            ind_test = opt.individual_test_specification
            sens = ind_test.sensitivity
            spec = ind_test.specificity
            test_lag = ind_test.test_result_lag

            alertthreshold = opt.alert_threshold
            accuracy = opt.accuracy

            return percent_clinic_tested, sens,
                spec, test_lag, alertthreshold, accuracy
        end
        reduce(vcat, _)
        DataFrames.DataFrame(
            _,
            [
                "percent_clinic_tested",
                "sensitivity",
                "specificity",
                "test_lag",
                "alert_threshold",
                "accuracy",
            ],
        )
    end
end

"""
    create_wide_optimal_thresholds_df(df, characteristic_sym)

Create a wide-format DataFrame for a specific characteristic.
"""
function create_wide_optimal_thresholds_df(df, characteristic_sym)
    maindf = Chain.@chain df begin
        DataFrames.select(
            _,
            DataFrames.Cols(
                x -> startswith(x, "s"),
                :test_lag,
                x -> contains(x, "tested"),
                characteristic_sym,
            ),
        )
        DataFrames.unstack(
            _,
            [:sensitivity, :specificity, :test_lag],
            :percent_clinic_tested,
            characteristic_sym,
        )
        DataFrames.select(
            _,
            DataFrames.Cols(
                x -> startswith(x, "s"),
                "test_lag",
                :,
            ),
        )
    end

    clinical_case_df = DataFrames.subset(
        maindf,
        :sensitivity => sens -> sens .== 1.0,
        :specificity => spec -> spec .== 0 .|| spec .== 0.8,
    )

    return Chain.@chain maindf begin
        DataFrames.antijoin(
            _, clinical_case_df; on = [:test_lag, :sensitivity, :specificity]
        )
        DataFramesMeta.@orderby :specificity, :test_lag
        vcat(clinical_case_df, _)
    end
end

"""
    create_optimal_threshold_summary_df(optimal_thresholds_vec, characteristic; 
                                       percentiles=[0.25, 0.5, 0.75], 
                                       nboots=10000, ci=0.95)

Create summary DataFrame for optimal threshold characteristics.
"""
function create_optimal_threshold_summary_df(
        optimal_thresholds_vec,
        characteristic::Symbol;
        percentiles = [0.25, 0.5, 0.75],
        nboots = 10000,
        ci = 0.95,
    )
    return create_optimal_threshold_summary_df(
        optimal_thresholds_vec,
        [characteristic];
        percentiles = percentiles,
        nboots = nboots,
        ci = ci,
    )
end

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
                    getfield.(opt.outbreak_threshold_chars, characteristic);
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
        "percent_clinic_tested" =>
            getfield.(optimal_thresholds_vec, :percent_clinic_tested),
        "sensitivity" =>
            getfield.(
            optimal_thresholds_vec.individual_test_specification,
            :sensitivity,
        ),
        "specificity" =>
            getfield.(
            optimal_thresholds_vec.individual_test_specification,
            :specificity,
        ),
        "test_lag" =>
            getfield.(
            optimal_thresholds_vec.individual_test_specification,
            :test_result_lag,
        ),
        "alert_threshold" =>
            getfield.(
            optimal_thresholds_vec,
            :alert_threshold,
        ),
        "accuracy" => getfield.(
            optimal_thresholds_vec,
            :accuracy,
        ),
    )

    return hcat(core_df, chars_df)
end

"""
    create_all_wide_optimal_threshold_summary_dfs(df; summary_stats=...)

Create wide-format DataFrames for all summary statistics.
"""
function create_all_wide_optimal_threshold_summary_dfs(
        df; summary_stats = ["mean", "lower_CI_95", "upper_CI_95", 0.25, 0.5, 0.75]
    )
    summary_stats_labels = summary_stats
    summary_stats_symbols = copy(summary_stats)

    for (i, stat) in pairs(summary_stats_labels)
        if typeof(stat) <: AbstractString
            summary_stats_symbols[i] = Symbol(stat)
            continue
        end
        summary_stats_labels[i] = "$(Int64(stat * 100))th"
        summary_stats_symbols[i] = Symbol("perc_$(Int64(stat * 100))th")
    end

    summary_dfs = map(
        char -> create_wide_optimal_threshold_summary_df(df, char),
        summary_stats_labels,
    )

    return (; zip(summary_stats_symbols, summary_dfs)...)
end

"""
    create_wide_optimal_threshold_summary_df(df, characteristic)

Create wide-format DataFrame for a specific summary characteristic.
"""
function create_wide_optimal_threshold_summary_df(df, characteristic)
    characteristic_str = filter(
        col -> contains(col, characteristic), names(df)
    )[1]
    characteristic_sym = Symbol(characteristic_str)

    return create_wide_optimal_thresholds_df(df, characteristic_sym)
end
