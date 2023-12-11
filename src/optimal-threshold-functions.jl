using StructArrays
using XLSX: XLSX

function calculate_OptimalThresholdCharacteristics(
    percent_clinic_tested_vec,
    ind_test_spec_vec,
    base_parameters
)
    non_clinical_case_test_spec_vec = filter(
        spec -> !(spec in CLINICAL_TEST_SPECS),
        ind_test_spec_vec
    )

    non_clinical_case_optimal_thresholds_vec = Vector{
        OptimalThresholdCharacteristics
    }(
        undef,
        length(percent_clinic_tested_vec) *
        length(non_clinical_case_test_spec_vec),
    )

    @showprogress for (i, (percent_clinic_tested, ind_test_spec)) in enumerate(
        Iterators.product(
            percent_clinic_tested_vec, non_clinical_case_test_spec_vec
        ),
    )
        non_clinical_case_optimal_thresholds_vec[i] = calculate_optimal_threshold(
            percent_clinic_tested,
            ind_test_spec,
            base_parameters
        )
    end

    clinical_case_optimal_thresholds_vec = Vector{
        OptimalThresholdCharacteristics
    }(
        undef, length(CLINICAL_TEST_SPECS)
    )

    for (i, ind_test_spec) in enumerate(CLINICAL_TEST_SPECS)
        clinical_case_optimal_thresholds_vec[i] = calculate_optimal_threshold(
            1.0,
            ind_test_spec,
            base_parameters
        )
    end

    return StructArray(
        vcat(
            non_clinical_case_optimal_thresholds_vec,
            clinical_case_optimal_thresholds_vec,
        ),
    )
end

function calculate_optimal_threshold(
    percent_clinic_tested,
    individual_test_specification,
    base_parameters
)
    @unpack alertthreshold_vec,
    ensemble_specification,
    noise_specification,
    outbreak_specification,
    moving_avg_detection_lag,
    percent_visit_clinic = base_parameters

    ensemble_scenario_spec_vec = map(
        threshold -> ScenarioSpecification(
            ensemble_specification,
            outbreak_specification,
            noise_specification,
            OutbreakDetectionSpecification(
                threshold,
                moving_avg_detection_lag,
                percent_visit_clinic,
                percent_clinic_tested,
            ),
            individual_test_specification,
        ),
        alertthreshold_vec,
    )

    optimal_accuracy = 0.0
    optimal_threshold = 0
    optimal_OT_chars = 0

    for (i, ensemble_scenario_spec) in pairs(ensemble_scenario_spec_vec)
        OT_chars = get_ensemble_file(ensemble_scenario_spec)["OT_chars"]
        accuracy = median(OT_chars.accuracy)
        alert_threshold =
            ensemble_scenario_spec.outbreak_detection_specification.alert_threshold

        if i == 1
            optimal_accuracy = accuracy
            optimal_threshold = alert_threshold
            optimal_OT_chars = OT_chars
            continue
        end
        if !isnan(accuracy) && accuracy > optimal_accuracy
            optimal_accuracy = accuracy
            optimal_threshold = alert_threshold
            optimal_OT_chars = OT_chars
        end
    end

    return OptimalThresholdCharacteristics(
        optimal_OT_chars,
        individual_test_specification,
        percent_clinic_tested,
        optimal_threshold,
        optimal_accuracy,
    )
end

function create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec,
    characteristic;
    percentiles = [0.25, 0.5, 0.75],
    filepath = datadir("optimal-threshold-results"),
    kwargs...,
)
    kwargs_dict = Dict(kwargs)

    long_df = create_optimal_threshold_summary_df(
        optimal_thresholds_vec, characteristic; percentiles = percentiles
    )

    if haskey(kwargs_dict, :cfrs)
        base_df = long_df[
            :,
            [
                "percent_clinic_tested",
                "sensitivity",
                "specificity",
                "test_lag",
                "alert_threshold",
                "accuracy",
            ],
        ]

        for (country, cfr) in pairs(kwargs_dict[:cfrs])
            cfr_df = long_df[:, Not(names(base_df))]
            cfr_df .*= cfr

            cfr_long_df = hcat(base_df, cfr_df)

            cfr_wide_df_tuples = create_all_wide_optimal_threshold_summary_dfs(
                cfr_long_df
            )

            cfr_info = DataFrame(; Country = String(country), CFR = cfr)

            round_cfr = round(cfr; digits = 5)

            cfr_filename = "optimal-threshold-result-tables_$(characteristic)_$(country)_CFR_$(round_cfr)"
            save_xlsx_optimal_threshold_summaries(
                (; cfr_info, cfr_long_df, cfr_wide_df_tuples...), cfr_filename;
                filepath = filepath,
            )
        end
    end

    wide_df_tuples = create_all_wide_optimal_threshold_summary_dfs(long_df)

    filename = "optimal-threshold-result-tables_$(characteristic)"
    save_xlsx_optimal_threshold_summaries(
        (; long_df, wide_df_tuples...), filename; filepath = filepath
    )

    @info "Saved the summary statistics for $(characteristic)"

    return nothing
end

function create_all_wide_optimal_threshold_summary_dfs(
    df; summary_stats = ["mean", 0.25, 0.50, 0.75]
)
    summary_stats_labels = summary_stats
    summary_stats_symbols = copy(summary_stats)

    for (i, stat) in pairs(summary_stats_labels)
        if stat == "mean"
            summary_stats_symbols[i] = Symbol(stat)
            continue
        end
        summary_stats_labels[i] = "$(Int64(stat*100))th"
        summary_stats_symbols[i] = Symbol("perc_$(Int64(stat*100))th")
    end

    summary_dfs = map(
        char -> create_wide_optimal_threshold_summary_df(df, char),
        summary_stats_labels,
    )

    return (; zip(summary_stats_symbols, summary_dfs)...)
end

function create_wide_optimal_threshold_summary_df(df, characteristic)
    characteristic_str = filter(
        col -> contains(col, characteristic), names(df)
    )[1]
    characteristic_sym = Symbol(characteristic_str)

    return create_wide_optimal_thresholds_df(df, characteristic_sym)
end

function create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec;
    filepath = datadir("optimal-threshold-results")
)
    long_df = create_optimal_thresholds_df(
        optimal_thresholds_vec
    )

    alert_thresholds = create_wide_optimal_thresholds_df(
        long_df, :alert_threshold
    )
    accuracy = create_wide_optimal_thresholds_df(
        long_df, :accuracy
    )

    filename = "optimal-threshold-result-tables_thresholds"
    save_xlsx_optimal_threshold_summaries(
        (; long_df, alert_thresholds, accuracy), filename; filepath = filepath
    )

    @info "Saved the thresholds and accuracy table"

    return nothing
end

function create_optimal_thresholds_df(optimal_thresholds_vec)
    @chain begin
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
        DataFrame(
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

function create_wide_optimal_thresholds_df(df, characteristic_sym)
    maindf = @chain df begin
        select(
            _,
            Cols(
                x -> startswith(x, "s"),
                :test_lag,
                x -> contains(x, "tested"),
                characteristic_sym,
            ),
        )
        unstack(
            _,
            [:sensitivity, :specificity, :test_lag],
            :percent_clinic_tested,
            characteristic_sym,
        )
        select(
            _,
            Cols(
                x -> startswith(x, "s"),
                "test_lag",
                x -> startswith(x, "0"),
                "1.0",
            ),
        )
    end

    clinical_case_df = subset(
        maindf,
        :sensitivity => sens -> sens .== 1.0,
        :specificity => spec -> spec .== 0 .|| spec .== 0.8,
    )

    return @chain maindf begin
        antijoin(
            _, clinical_case_df; on = [:test_lag, :sensitivity, :specificity]
        )
        @orderby :specificity, :test_lag
        vcat(clinical_case_df, _)
    end
end

function create_optimal_threshold_summary_df(
    optimal_thresholds_vec,
    characteristic;
    percentiles = [0.25, 0.5, 0.75]
)
    percentile_labels = map(perc -> "$(Int64(perc*100))th", percentiles)
    char_percentile_labels = map(
        perc -> "$(characteristic)_$perc", percentile_labels
    )
    char_mean_label = "$(characteristic)_mean"

    @chain begin
        map(optimal_thresholds_vec) do opt
            percent_clinic_tested = opt.percent_clinic_tested

            ind_test = opt.individual_test_specification
            sens = ind_test.sensitivity
            spec = ind_test.specificity
            test_lag = ind_test.test_result_lag

            alertthreshold = opt.alert_threshold
            accuracy = opt.accuracy

            char_mean, chars_percentiles = calculate_optimal_threshold_summaries(
                getfield.(opt.outbreak_threshold_chars, characteristic);
                percentiles = percentiles,
            )

            return percent_clinic_tested, sens,
            spec, test_lag, alertthreshold, accuracy,
            char_mean, chars_percentiles...
        end
        reduce(vcat, _)
        DataFrame(
            _,
            [
                "percent_clinic_tested",
                "sensitivity",
                "specificity",
                "test_lag",
                "alert_threshold",
                "accuracy",
                char_mean_label,
                char_percentile_labels...,
            ],
        )
    end
end

function calculate_optimal_threshold_summaries(
    char_vecs; percentiles = [0.25, 0.5, 0.75]
)
    all_chars = reduce(vcat, char_vecs)

    char_percentiles = map(
        percentile -> quantile(all_chars, percentile), percentiles
    )
    char_mean = mean(all_chars)

    return char_mean, char_percentiles
end

function save_xlsx_optimal_threshold_summaries(
    summary_tuple, filename; filepath = datadir("optimal-threshold-results")
)
    sheet_names = String.(keys(summary_tuple))
    file = joinpath(filepath, "$(filename).xlsx")

    XLSX.openxlsx(file; mode = "w") do xf
        for i in eachindex(sheet_names)
            sheet_name = sheet_names[i]
            df = summary_tuple[i]

            if i == firstindex(sheet_names)
                sheet = xf[1]
                XLSX.rename!(sheet, sheet_name)
                XLSX.writetable!(sheet, df)
            else
                sheet = XLSX.addsheet!(xf, sheet_name)
                XLSX.writetable!(sheet, df)
            end
        end
    end
end
