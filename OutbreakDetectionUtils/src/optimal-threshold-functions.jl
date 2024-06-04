using ProgressMeter: ProgressMeter
using StructArrays: StructArrays
using DataFrames: DataFrames
using DataFramesMeta: DataFramesMeta
using Chain: Chain
using UnPack: UnPack
using Match: Match
using XLSX: XLSX
using RCall: RCall

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
        OptimalThresholdCharacteristics
    }(
        undef,
        length(percent_clinic_tested_vec) *
        length(non_clinical_case_test_spec_vec),
    )

    ProgressMeter.@showprogress for (
        i, (percent_clinic_tested, ind_test_spec)
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

    clinical_case_optimal_thresholds_vec = Vector{
        OptimalThresholdCharacteristics
    }(
        undef, length(CLINICAL_TEST_SPECS)
    )

    for (i, ind_test_spec) in enumerate(CLINICAL_TEST_SPECS)
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
    tabledirpath = outdir("ensemble/optimal-threshold-results"),
    filename = "",
    kwargs...,
)
    kwargs_dict = Dict(kwargs)

    long_df = create_optimal_threshold_summary_df(
        optimal_thresholds_vec, characteristic; percentiles = percentiles
    )

    base_columns = [
        "percent_clinic_tested",
        "sensitivity",
        "specificity",
        "test_lag",
        "alert_threshold",
        "accuracy",
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
                country_long_df
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

            if haskey(kwargs_dict, :gt_kwargs)
                gt_kwargs = kwargs_dict[:gt_kwargs]

                if !haskey(gt_kwargs, :summary_stats)
                    @error "gt_kwargs does not have summary_stats. Please provide one of the following: mean, perc_25th, perc_50th, perc_75th (or any other percentile that was calculated)"
                end

                map(gt_kwargs.summary_stats) do stat
                    statdf = getproperty(country_wide_df_tuples, Symbol(stat))

                    gt_table(
                        statdf;
                        testing_rates = gt_kwargs.testing_rates,
                        colorschemes = gt_kwargs.colorschemes,
                        save = gt_kwargs.save,
                        show = gt_kwargs.show,
                        filepath = tabledirpath,
                        filename = country_filename * "_$(stat).png",
                        decimals = gt_kwargs.decimals,
                        gt_kwargs...,
                    )
                end
            end

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
                    cfr_long_df
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
        wide_df_tuples = create_all_wide_optimal_threshold_summary_dfs(long_df)

        if haskey(kwargs_dict, :gt_kwargs)
            gt_kwargs = kwargs_dict[:gt_kwargs]

            if !haskey(gt_kwargs, :summary_stats)
                @error "gt_kwargs does not have summary_stats. Please provide one of the following: mean, perc_25th, perc_50th, perc_75th (or any other percentile that was calculated)"
            end

            map(gt_kwargs.summary_stats) do stat
                statdf = getproperty(wide_df_tuples, Symbol(stat))

                gt_table(
                    statdf;
                    testing_rates = gt_kwargs.testing_rates,
                    colorschemes = gt_kwargs.colorschemes,
                    save = gt_kwargs.save,
                    show = gt_kwargs.show,
                    filepath = tabledirpath,
                    filename = base_filename * "_$(stat).png",
                    decimals = gt_kwargs.decimals,
                    gt_kwargs...,
                )
            end
        end

        save_xlsx_optimal_threshold_summaries(
            (; long_df, wide_df_tuples...), base_filename;
            filepath = tabledirpath,
        )
    end

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

    if haskey(kwargs_dict, :gt_kwargs)
        gt_kwargs = kwargs_dict[:gt_kwargs]

        alert_thresholds_kwargs = (;
            testing_rates = gt_kwargs.testing_rates,
            colorschemes = gt_kwargs.alert_threshold_colorscheme,
            filepath = tabledirpath,
            filename = filename * "_alert_thresholds.png",
            save = gt_kwargs.save,
            show = gt_kwargs.show,
            decimals = 0,
        )

        if haskey(gt_kwargs, :alert_threshold_domain)
            alert_thresholds_kwargs = (;
                alert_thresholds_kwargs...,
                domain = gt_kwargs.alert_threshold_domain,
            )
        end

        gt_table(alert_thresholds; alert_thresholds_kwargs...)

        accuracy_kwargs = (;
            testing_rates = gt_kwargs.testing_rates,
            colorschemes = gt_kwargs.accuracy_colorscheme,
            filepath = tabledirpath,
            filename = filename * "_accuracy.png",
            save = gt_kwargs.save,
            show = gt_kwargs.show,
        )

        if haskey(gt_kwargs, :accuracy_domain)
            accuracy_kwargs = (;
                accuracy_kwargs...,
                domain = gt_kwargs.accuracy_domain,
            )
        end

        gt_table(accuracy; accuracy_kwargs...)
    end

    @info "Saved the thresholds and accuracy table"

    return nothing
end

function create_optimal_thresholds_df(optimal_thresholds_vec)
    Chain.@chain begin
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
    characteristic;
    percentiles = [0.25, 0.5, 0.75],
)
    percentile_labels = map(perc -> "$(Int64(perc*100))th", percentiles)
    char_percentile_labels = map(
        perc -> "$(characteristic)_$perc", percentile_labels
    )
    char_mean_label = "$(characteristic)_mean"

    Chain.@chain begin
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
        DataFrames.DataFrame(
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

    if isempty(all_chars)
        return missing, repeat([missing], length(percentiles))
    end

    char_percentiles = map(
        percentile -> quantile(all_chars, percentile), percentiles
    )
    char_mean = mean(all_chars)

    return char_mean, char_percentiles
end

function save_xlsx_optimal_threshold_summaries(
    summary_tuple, filename; filepath = outdir("optimal-threshold-results")
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

function _rename_test_scenario(x)
    Match.@match x begin
        0.85 => "RDT Equivalent (0.85)"
        0.9 => "RDT Equivalent (0.9)"
        _ => "ELISA Equivalent"
    end
end

function gt_table(
    df;
    testing_rates = DataFrames.Between("0.1", "0.6"),
    colorschemes = ["ggsci::green_material"],
    save = "no",
    show = "yes",
    filepath = outdir("tables"),
    filename = "optimal_thresholds.png",
    decimals = 2,
    container_width_px = 960,
    container_height_px = 660,
    table_width_pct = 100,
    kwargs...,
)
    kwarg_dict = Dict(kwargs...)

    filtered = Chain.@chain df begin
        DataFramesMeta.@rsubset(:specificity > 0.8)
        DataFramesMeta.@rtransform(
            :test_scenario = _rename_test_scenario(:sensitivity)
        )
        DataFrames.select(:test_scenario, :test_lag, testing_rates)
        DataFrames.rename(:test_scenario => "Test Scenario", :test_lag => "Lag")
    end

    matrix = Matrix(filtered[:, testing_rates])
    maxval = maximum(matrix)
    minval = minimum(matrix)

    if !haskey(kwarg_dict, :domain)
        domain = (minval, maxval)
    else
        domain = kwarg_dict[:domain]
    end

    if length(colorschemes) !== 1 && !haskey(kwarg_dict, :domain)
        if length(colorschemes) > 2
            @error "More than 2 colorschemes provided"
            @show colorschemes
            @show length(colorschemes)
        end
        domain = ((minval, 0), (0, maxval))
    end

    mkpath(filepath)

    RCall.@rput filtered save show filepath filename domain colorschemes decimals container_width_px container_height_px table_width_pct

    RCall.R"""
    library(gt)
    library(tidyverse)

    table <- filtered %>%
     gt() %>%
     tab_spanner(label = "Test Characteristic", columns = 1:2) %>%
     tab_spanner(label = "Testing Rate", columns = 3:ncol(filtered)) %>%
     fmt_number(columns = 3:ncol(filtered), decimals = decimals) %>%
     opt_table_font(
        font = list(
            google_font("Lato"),
            default_fonts()
        ),
        weight = 300
     ) %>%
     tab_options(
        table.font.size = gt::px(23L),
        table.width = pct(table_width_pct),
        container.width = px(container_width_px),
        container.height = px(container_height_px),
     ) %>%
     tab_style(
        style = cell_text(weight = 900),
        locations = list(
            cells_column_spanners(spanners = everything()),
            cells_column_labels(columns = everything())
        )
     ) %>%
     cols_width(
        2:ncol(filtered) ~ pct(65/(ncol(filtered) - 2))
    ) %>%
    cols_align(align = "center", columns = 2:ncol(filtered))


    if (length(colorschemes) == 1) {
    table <- table %>%
        data_color(
            columns = 3:ncol(filtered),
            domain = domain,
            palette = colorschemes
        )
    } else {
        colorpalette <- function(x) {
          f_neg <- scales::col_numeric(
            palette = c(paletteer::paletteer_d(colorschemes[1])[9], '#ffffff'),
            domain = domain[1],
          )
          f_pos <- scales::col_numeric(
            palette = c('#ffffff', paletteer::paletteer_d(colorschemes[2])[9]),
            domain = domain[2]
          )
          ifelse(x < 0, f_neg(x), f_pos(x))
        }

        table <- table %>%
            data_color(
                columns = 3:ncol(filtered),
                fn = colorpalette
            )
    }

    if (save == "yes") {
        gt::gtsave(table, filename, filepath)
    }

    if (show == "yes") {
        table
    }

    if (show == "no" && save == "no") {
        print("The table should be saved or shown.")
    }
    """
end
