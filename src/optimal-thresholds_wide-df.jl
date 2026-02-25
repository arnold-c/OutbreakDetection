export create_wide_optimal_thresholds_df

"""
    create_wide_optimal_thresholds_df(optimal_thresholds_vec, characteristic_sym;
                                      simplify=false)

Create a wide table directly from optimization results.
"""
function create_wide_optimal_thresholds_df(
        optimal_thresholds_vec::StructVector{OutbreakDetectionCore.OptimizationResult},
        characteristic_sym::Symbol;
        simplify::Bool = false,
    )
    possible_properties = filter(
        p -> !(
            p in [
                :ensemble_specification,
                :noise_level,
                :noise_type_description,
                :vaccination_coverage,
                :test_specification,
                :percent_tested,
                :alert_method,
                :accuracy_metric,
                :threshold_bounds,
                :outbreak_specification,
                :alert_filtering_strategy,
            ]
        ),
        propertynames(optimal_thresholds_vec[1])
    )

    @assert characteristic_sym in possible_properties "characteristic_sym kwarg :$(characteristic_sym) is not a property name of struct OptimizationResult. Possible options are $(possible_properties)"

    long_df = DataFrames.DataFrame(
        (
                noise_type_description = opt.noise_type_description,
                noise_level = opt.noise_level,
                percent_tested = opt.percent_tested,
                sensitivity = opt.test_specification.sensitivity,
                specificity = opt.test_specification.specificity,
                test_lag = opt.test_specification.test_result_lag,
                characteristic_value = _characteristic_value_for_table(
                    opt,
                    characteristic_sym,
                ),
            ) for opt in optimal_thresholds_vec
    )

    return create_wide_optimal_thresholds_df(
        long_df,
        :characteristic_value;
        simplify = simplify,
    )
end

"""
    create_wide_optimal_thresholds_df(df, characteristic_sym)

Create a wide-format DataFrame for a specific characteristic.
"""
function create_wide_optimal_thresholds_df(
        df::DataFrames.DataFrame,
        characteristic_col::Symbol;
        simplify::Bool = false,
    )

    DataFrames.rename!(
        df,
        :noise_type_description => :noise_type
    )

    id_columns = [
        :noise_type,
        :noise_level,
        :sensitivity,
        :specificity,
        :test_lag,
    ]

    working_df = DataFrames.select(
        df,
        unique(vcat(id_columns, [:percent_tested, characteristic_col])),
    )

    wide_df = DataFrames.unstack(
        working_df,
        id_columns,
        :percent_tested,
        :characteristic_value,
    )

    wide_df[!, :test_type] = table_test_type.(
        Float64.(wide_df.sensitivity),
        Float64.(wide_df.specificity),
        Int64.(wide_df.test_lag),
    )

    sort!(
        wide_df,
        [
            :noise_type,
            :noise_level,
            :sensitivity,
            :specificity,
            DataFrames.order(:test_lag; rev = true),
        ]
    )

    percent_columns = filter(_is_percent_tested_column, propertynames(wide_df))
    sort!(percent_columns; by = x -> parse(Float64, string(x)))

    leading_columns = filter(
        col -> col in propertynames(wide_df),
        [:noise_type, :noise_level, :test_type, :test_lag],
    )

    output_df = DataFrames.select(
        wide_df,
        vcat(leading_columns, percent_columns)
    )

    if simplify
        return _simplify_wide_optimal_thresholds_df(output_df)
    end

    return output_df
end

"""
    _simplify_wide_optimal_thresholds_df(wide_df)

Collapse duplicated perfect-test rows across noise levels and enforce table order:
static, dynamic, then all-noise perfect rows.
"""
function _simplify_wide_optimal_thresholds_df(wide_df::DataFrames.DataFrame)
    required_columns = (:test_type, :noise_type, :noise_level)

    percent_columns = filter(_is_percent_tested_column, propertynames(wide_df))

    perfect_mask = wide_df.test_type .== "Perfect Test"

    non_perfect_df = wide_df[.!perfect_mask, :]
    perfect_df = wide_df[perfect_mask, :]

    perfect_collapsed_df = DataFrames.combine(
        DataFrames.groupby(perfect_df, [:test_type, :test_lag]),
        [column_name => first => column_name for column_name in percent_columns],
    )

    perfect_collapsed_df[!, :noise_type] .= :all_noise_structures
    perfect_collapsed_df[!, :noise_level] .= missing

    DataFrames.select!(perfect_collapsed_df, propertynames(wide_df))

    merged_df = vcat(non_perfect_df, perfect_collapsed_df; cols = :union)

    merged_df[!, :noise_sort_order] = map(merged_df.noise_type) do noise_type
        if noise_type == :static
            return 1
        elseif noise_type == :dynamic
            return 2
        elseif noise_type == :all_noise_structures
            return 3
        end

        return 99
    end

    sort!(merged_df, [:noise_sort_order])
    DataFrames.select!(merged_df, DataFrames.Not(:noise_sort_order))

    return merged_df
end
