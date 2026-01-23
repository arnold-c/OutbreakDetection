"""
    reshape_optimization_results_to_matrix(results::StructVector{OptimizationResult})

Reshape a StructVector of OptimizationResult into a matrix organized by noise type and magnitude.
Returns a matrix where rows represent noise types (static vs dynamical) and columns represent
noise magnitudes/scalings.
"""
function reshape_optimization_results_to_matrix(
        results::StructVector{OutbreakDetectionCore.OptimizationResult}
    )
    # Get unique noise types and levels
    unique_noise_types = unique(results.noise_type_description)

    # Group results by noise type
    noise_type_groups = map(unique_noise_types) do noise_type
        filter(r -> r.noise_type_description == noise_type, results)
    end

    # For each noise type group, sort by noise level and get unique levels
    sorted_groups = map(noise_type_groups) do group
        unique_levels = sort(unique(group.noise_level))
        return (group = group, levels = unique_levels)
    end

    # Determine matrix dimensions
    num_noise_types = length(unique_noise_types)
    num_noise_levels = length(sorted_groups[1].levels)

    # Create matrix to hold StructVectors for each cell
    result_matrix = Matrix{StructVector{OutbreakDetectionCore.OptimizationResult}}(
        undef, num_noise_types, num_noise_levels
    )

    # Fill matrix
    for (i, (noise_type, sorted_group)) in enumerate(zip(unique_noise_types, sorted_groups))
        for (j, noise_level) in enumerate(sorted_group.levels)
            # Filter results for this specific noise type and level
            cell_results = filter(
                r -> r.noise_type_description == noise_type && r.noise_level == noise_level,
                sorted_group.group
            )

            # Sort by test specification for consistent ordering
            # Sort by: test_result_lag (descending), then sensitivity (ascending), then specificity (ascending)
            test_specs = cell_results.test_specification
            sort_order = sortperm(
                collect(enumerate(test_specs));
                by = x -> (-x[2].test_result_lag, x[2].sensitivity, x[2].specificity)
            )
            sort_indices = [x[1] for x in sort_order]

            result_matrix[i, j] = cell_results[sort_indices]
        end
    end

    return result_matrix, unique_noise_types
end

"""
    compute_summary_statistics(values::Vector{T}; percentiles=[0.1, 0.9]) where T

Compute mean and percentiles for a vector of values.
Returns a NamedTuple with mean and percentile values.
"""
function compute_summary_statistics(
        values::Vector{T};
        percentiles = [0.1, 0.9]
    ) where {T}
    if isempty(values)
        return (mean = NaN, percentiles = fill(NaN, length(percentiles)))
    end

    mean_val = StatsBase.mean(values)
    percentile_vals = StatsBase.percentile(values, percentiles .* 100)

    return (mean = mean_val, percentiles = percentile_vals)
end

"""
    compute_nested_summary_statistics(nested_values::Vector{Vector{T}}; percentiles=[0.1, 0.9]) where T

Compute mean and percentiles for nested vectors (e.g., detection_delays).
First flattens the nested structure, then computes statistics.
"""
function compute_nested_summary_statistics(
        nested_values::Vector{Vector{T}};
        percentiles = [0.1, 0.9]
    ) where {T}
    if isempty(nested_values)
        return (mean = NaN, percentiles = fill(NaN, length(percentiles)))
    end

    # Flatten nested vectors
    flattened = reduce(vcat, nested_values)

    return compute_summary_statistics(flattened; percentiles = percentiles)
end

"""
    extract_outcome_values(result::OptimizationResult, outcome::Symbol)

Extract the appropriate field from OptimizationResult based on the outcome symbol.
"""
function extract_outcome_values(
        result::OutbreakDetectionCore.OptimizationResult,
        outcome::Symbol
    )
    if outcome == :accuracy
        return result.accuracies
    elseif outcome == :alert_threshold
        return result.optimal_threshold
    elseif outcome == :detectiondelays
        return result.detection_delays
    elseif outcome == :proportion_timeseries_in_alert
        return result.proportion_timeseries_in_alert
    elseif outcome == :proportion_timeseries_in_outbreak
        return result.proportion_timeseries_in_outbreak
    elseif outcome == :proportion_alerts_correct
        return result.proportion_alerts_correct
    elseif outcome == :proportion_outbreaks_detected
        return result.proportion_outbreaks_detected
    elseif outcome == :unavoidable_cases
        return result.unavoidable_cases
    elseif outcome == :avoidable_cases
        # Compute avoidable cases from unavoidable cases
        # This would need the total cases, which we don't have directly
        # For now, return unavoidable_cases as placeholder
        @warn "avoidable_cases computation not yet implemented"
        return result.unavoidable_cases
    elseif outcome == :alert_durations
        return result.alert_durations
    elseif outcome == :outbreak_durations
        return result.outbreak_durations
    else
        error("Unknown outcome: $outcome")
    end
end

"""
    get_noise_label(noise_type::Symbol)

Convert noise type symbol to human-readable label.
"""
function get_noise_label(noise_type::Symbol)
    if noise_type == :poisson || noise_type == :static
        return "Static Noise"
    elseif noise_type == :dynamical_inphase || noise_type == :dynamic
        return "Dynamical Noise"
    else
        return titlecase(string(noise_type))
    end
end
