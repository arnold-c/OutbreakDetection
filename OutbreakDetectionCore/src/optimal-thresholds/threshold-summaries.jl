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
