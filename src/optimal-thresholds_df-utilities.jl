export prepare_wide_optimal_thresholds_df_format,
    subset_for_noise_level_with_perfect_tests

"""
    prepare_wide_optimal_thresholds_df_format(wide_thresholds_df)

Prepare a wide optimal-threshold table for display by:
- converting `noise_type` labels to human-readable strings,
- renaming identifier columns for publication-ready headers, and
- converting numeric percent-tested column names (for example, `0.2`) into
  percentage labels (for example, `"20.0%"`).
"""
function prepare_wide_optimal_thresholds_df_format(wide_thresholds_df)
    df = copy(wide_thresholds_df)
    df.noise_type .= String.(df.noise_type)
    replace!(
        df.noise_type,
        "static" => "Static Noise",
        "dynamic" => "Dynamic Noise",
        "all_noise_structures" => "All noise structures",
    )

    DataFrames.rename!(
        df,
        :noise_type => "Noise Type",
        :test_type => "Test Type",
        :test_lag => "Test Lag",
    )

    DataFrames.rename!(
        prop -> Printf.@sprintf("%i%%", parse(Float64, prop) * 100),
        df;
        cols = names(df, r"\d")
    )
    return df
end

"""
    subset_for_noise_level_with_perfect_tests(wide_df, noise_level)

Subset a wide optimal-threshold table to a single `noise_level` while keeping the
collapsed perfect-test rows (`noise_type == :all_noise_structures`).

This is intended for tables generated with
`create_wide_optimal_thresholds_df(...; simplify=true)`.
"""
function subset_for_noise_level_with_perfect_tests(wide_df, noise_level)
    return DataFrames.subset(
        wide_df,
        [:noise_level, :noise_type] => (levels, noise_types) ->
        ((levels .== noise_level) .| (noise_types .== :all_noise_structures)),
    )
end

"""
    _is_percent_tested_column(column_name)

Return true if `column_name` parses to a numeric testing proportion.
"""
function _is_percent_tested_column(column_name)
    try
        parse(Float64, string(column_name))
        return true
    catch
        return false
    end
end

"""
    _characteristic_value_for_table(opt, characteristic_sym)

Return a string-formatted characteristic value for table output.

If the selected characteristic is a nested vector, values are flattened and
averaged before formatting to two decimal places.
"""
function _characteristic_value_for_table(opt, characteristic_sym::Symbol)
    characteristic_value = getproperty(opt, characteristic_sym)

    result = characteristic_value
    if characteristic_value isa AbstractVector
        result = StatsBase.mean(Iterators.flatten(characteristic_value))
    end

    return @sprintf("%.2f", result)
end
