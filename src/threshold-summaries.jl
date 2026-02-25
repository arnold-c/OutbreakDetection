export calculate_optimal_threshold_summaries

"""
    calculate_optimal_threshold_summaries(char_vecs; percentiles=[0.25, 0.5, 0.75],
                                          nboots=10000, ci=0.95)

Calculate summary statistics for optimal threshold characteristics.

Returns lower CI, mean, upper CI, and percentiles.
"""
function calculate_optimal_threshold_summaries(
        char_vecs; percentiles = [0.25, 0.5, 0.75], nboots = 10000, ci = 0.95
    )
    all_chars = _flatten_summary_values(char_vecs)

    if isempty(all_chars)
        missing_percentiles = isnothing(percentiles) ? Any[] : fill(missing, length(percentiles))
        return missing, missing, missing, missing_percentiles
    end

    if !isnothing(percentiles)
        char_percentiles = map(
            percentile -> StatsBase.quantile(all_chars, percentile), percentiles
        )
    else
        char_percentiles = Any[]
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

"""
    _flatten_summary_values(values)

Flatten scalar/vector/nested-vector values into a single vector.
"""
function _flatten_summary_values(values)
    if values isa AbstractVector{<:AbstractVector}
        if isempty(values)
            return Any[]
        end
        return reduce(vcat, values)
    end

    if values isa AbstractVector
        return collect(values)
    end

    return [values]
end
