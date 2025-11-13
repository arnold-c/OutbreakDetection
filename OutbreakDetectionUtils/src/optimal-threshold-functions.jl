using ProgressMeter: ProgressMeter
using StructArrays: StructArrays
using DataFrames: DataFrames
using DataFramesMeta: DataFramesMeta
using Chain: Chain
using UnPack: UnPack
using StatsBase: StatsBase
using Bootstrap: Bootstrap
using Match: Match
using XLSX: XLSX

function calculate_OptimalThresholdCharacteristics(
        percent_clinic_tested_vec,
        ind_test_spec_vec,
        base_parameters,
    )
    non_clinical_case_test_spec_vec = filter(
        spec -> !(spec in CLINICAL_TEST_SPECS),
        ind_test_spec_vec,
    )

    non_clinical_case_optimal_thresholds_vec = Vector{
        OptimalThresholdCharacteristics,
    }(
        undef,
        length(percent_clinic_tested_vec) *
            length(non_clinical_case_test_spec_vec),
    )

    ProgressMeter.@showprogress for (
            i, (percent_clinic_tested, ind_test_spec),
        ) in enumerate(
            Iterators.product(
                percent_clinic_tested_vec, non_clinical_case_test_spec_vec
            ),
        )
        non_clinical_case_optimal_thresholds_vec[i] = calculate_optimal_threshold(
            percent_clinic_tested,
            ind_test_spec,
            base_parameters,
        )
    end

    clinical_case_test_spec_vec = filter(
        spec -> spec in CLINICAL_TEST_SPECS,
        ind_test_spec_vec,
    )

    if length(clinical_case_test_spec_vec) == 0
        return StructArrays.StructArray(
            non_clinical_case_optimal_thresholds_vec
        )
    end

    clinical_case_optimal_thresholds_vec = Vector{
        OptimalThresholdCharacteristics,
    }(
        undef,
        length(clinical_case_test_spec_vec),
    )

    for (i, ind_test_spec) in enumerate(clinical_case_test_spec_vec)
        clinical_case_optimal_thresholds_vec[i] = calculate_optimal_threshold(
            1.0,
            ind_test_spec,
            base_parameters,
        )
    end

    return StructArrays.StructArray(
        vcat(
            non_clinical_case_optimal_thresholds_vec,
            clinical_case_optimal_thresholds_vec,
        ),
    )
end

function calculate_optimal_threshold(
        percent_clinic_tested,
        individual_test_specification,
        base_parameters,
    )
    UnPack.@unpack alertthreshold_vec,
        ensemble_specification,
        noise_specification,
        outbreak_specification,
        moving_avg_detection_lag,
        percent_visit_clinic,
        alertmethod = base_parameters

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
                alertmethod,
            ),
            individual_test_specification,
        ),
        alertthreshold_vec,
    )

    optimal_accuracy = 0.0
    optimal_threshold = 0
    optimal_OT_chars = 0

    for (i, ensemble_scenario_spec) in pairs(ensemble_scenario_spec_vec)
        scenario_chars_file = get_ensemble_file(ensemble_scenario_spec)
        OT_chars = scenario_chars_file["OT_chars"]
        accuracy = StatsBase.median(OT_chars.accuracy)
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
        noise_specification,
        percent_clinic_tested,
        optimal_threshold,
        optimal_accuracy,
    )
end

function create_and_save_xlsx_optimal_threshold_summaries(
        optimal_thresholds_vec,
        characteristic;
        percentiles = [0.25, 0.5, 0.75],
        nboots = 10000,
        ci = 0.95,
        tabledirpath = outdir("ensemble/optimal-threshold-results"),
        filename = "",
        kwargs...,
    )
    kwargs_dict = Dict(kwargs)

    long_df = create_optimal_threshold_summary_df(
        optimal_thresholds_vec,
        characteristic;
        percentiles = percentiles,
        nboots = nboots,
        ci = ci,
    )

    base_columns = [
        "percent_clinic_tested",
        "sensitivity",
        "specificity",
        "test_lag",
        "alert_threshold",
        "accuracy",
    ]

    summary_stats = [
        "mean",
        "lower_CI_$(Int64(ci * 100))",
        "upper_CI_$(Int64(ci * 100))",
        0.25,
        0.5,
        0.75,
    ]

    mkpath(tabledirpath)
    base_filename = "$(filename)_$(characteristic)"

    if haskey(kwargs_dict, :scale_annual)
        DataFrames.transform!(
            long_df,
            DataFrames.Not(base_columns) .=>
                x -> x .* kwargs_dict[:scale_annual];
            renamecols = false,
        )

        base_filename = base_filename * "_annual_scale"
    end

    if haskey(kwargs_dict, :countries)
        for country in kwargs_dict[:countries]
            if !haskey(country, :scale_population)
                @error "Country $(country) has no scale_population. Please provide one"
            end

            country_long_df = DataFrames.transform(
                long_df,
                DataFrames.Not(base_columns) .=>
                    x -> x .* country.scale_population;
                renamecols = false,
            )

            country_wide_df_tuples = create_all_wide_optimal_threshold_summary_dfs(
                country_long_df;
                summary_stats = summary_stats,
            )

            country_info_df = DataFrames.DataFrame(
                hcat(country...), [keys(country)...]
            )

            if !haskey(country, :code)
                @error "Country $(country) has no code. Please provide one"
            end

            if !haskey(country, :year)
                @error "Country $(country) has no year. Please provide one"
            end

            country_filename =
                base_filename * "_$(country.code)_$(country.year)"

            save_xlsx_optimal_threshold_summaries(
                (; country_info_df, country_long_df, country_wide_df_tuples...),
                country_filename;
                filepath = tabledirpath,
            )

            if haskey(country, :cfr) &&
                    occursin("case", String(characteristic))
                cfr_long_df = DataFrames.transform(
                    country_long_df,
                    DataFrames.Not(base_columns) .=> x -> x .* country.cfr;
                    renamecols = false,
                )

                DataFrames.rename!(
                    name -> replace(name, r"case" => "death"), cfr_long_df
                )

                cfr_wide_df_tuples = create_all_wide_optimal_threshold_summary_dfs(
                    cfr_long_df;
                    summary_stats = summary_stats,
                )

                round_cfr = round(country.cfr; digits = 3)

                cfr_filename = country_filename * "_CFR_$(round_cfr)"
                save_xlsx_optimal_threshold_summaries(
                    (; country_info_df, cfr_long_df, cfr_wide_df_tuples...),
                    cfr_filename;
                    filepath = tabledirpath,
                )
            end
        end
    else
        wide_df_tuples = create_all_wide_optimal_threshold_summary_dfs(
            long_df;
            summary_stats = summary_stats,
        )

        save_xlsx_optimal_threshold_summaries(
            (; long_df, wide_df_tuples...), base_filename;
            filepath = tabledirpath,
        )
    end

    @info "Saved the summary statistics for $(characteristic)"

    return nothing
end

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

function create_wide_optimal_threshold_summary_df(df, characteristic)
    characteristic_str = filter(
        col -> contains(col, characteristic), names(df)
    )[1]
    characteristic_sym = Symbol(characteristic_str)

    return create_wide_optimal_thresholds_df(df, characteristic_sym)
end

function create_and_save_xlsx_optimal_threshold_summaries(
        optimal_thresholds_vec;
        tabledirpath = outdir("optimal-threshold-results"),
        filename = "optimal-threshold-result-tables",
        kwargs...,
    )
    kwargs_dict = Dict(kwargs)

    long_df = create_optimal_thresholds_df(
        optimal_thresholds_vec
    )

    alert_thresholds = create_wide_optimal_thresholds_df(
        long_df, :alert_threshold
    )
    accuracy = create_wide_optimal_thresholds_df(
        long_df, :accuracy
    )

    mkpath(tabledirpath)

    save_xlsx_optimal_threshold_summaries(
        (; long_df, alert_thresholds, accuracy), filename;
        filepath = tabledirpath,
    )

    @info "Saved the thresholds and accuracy table"

    return nothing
end

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

function save_xlsx_optimal_threshold_summaries(
        summary_tuple, filename; filepath = outdir("optimal-threshold-results")
    )
    sheet_names = String.(keys(summary_tuple))
    file = joinpath(filepath, "$(filename).xlsx")

    return XLSX.openxlsx(file; mode = "w") do xf
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

function _rename_test_scenario(x)
    return Match.@match x begin
        0.85 => "RDT Equivalent (0.85)"
        0.9 => "RDT Equivalent (0.9)"
        _ => "ELISA Equivalent"
    end
end
