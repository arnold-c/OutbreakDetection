export subset_for_noise_level_with_perfect_tests

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
    _is_percent_column(column_name)

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

function _characteristic_value_for_table(opt, characteristic_sym::Symbol)
    characteristic_value = getproperty(opt, characteristic_sym)

    result = characteristic_value
    if characteristic_value isa AbstractVector{<:AbstractVector}
        result = StatsBase.mean(Iterators.flatten(characteristic_value))
    end

    return @sprintf("%.2f", result)
end
