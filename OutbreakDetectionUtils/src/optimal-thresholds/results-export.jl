using DataFrames: DataFrames
using XLSX: XLSX
using Match: Match

export save_xlsx_optimal_threshold_summaries,
    create_and_save_xlsx_optimal_threshold_summaries

"""
    save_xlsx_optimal_threshold_summaries(summary_tuple, filename; filepath=...)

Save optimal threshold summaries to an Excel file.
"""
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

"""
    create_and_save_xlsx_optimal_threshold_summaries(optimal_thresholds_vec; 
                                                     tabledirpath=..., 
                                                     filename=..., kwargs...)

Create and save optimal threshold summary tables to Excel.
"""
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

"""
    create_and_save_xlsx_optimal_threshold_summaries(optimal_thresholds_vec, 
                                                     characteristic; 
                                                     percentiles=..., 
                                                     nboots=..., ci=..., 
                                                     tabledirpath=..., 
                                                     filename=..., kwargs...)

Create and save detailed optimal threshold summary statistics with bootstrapped CIs.
"""
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

function _rename_test_scenario(x)
    return Match.@match x begin
        0.85 => "RDT Equivalent (0.85)"
        0.9 => "RDT Equivalent (0.9)"
        _ => "ELISA Equivalent"
    end
end
