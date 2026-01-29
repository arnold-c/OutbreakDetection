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
    unique_noise_types = unique(results.noise_type_description) |>
        # sort so static noise appears first for reshaped matrix
        noise_description -> sort(noise_description; rev = true)

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

            # Sort results by test specification and percent_tested for consistent ordering
            # Uses the same ordering as sort_test_specifications:
            # test_result_lag (descending), then (sensitivity, specificity) (ascending),
            # then percent_tested (ascending)
            sorted_indices = sort_structvector_by_test_specification_and_percent(
                cell_results.test_specification,
                cell_results.percent_tested
            )

            sorted_results = cell_results[sorted_indices]

            result_matrix[i, j] = sorted_results
        end
    end

    return result_matrix, unique_noise_types
end

"""
    sort_structvector_by_test_specification_and_percent(
        test_specifications,
        percent_tested
    )

Sort indices of optimization results by test specification and then by percent_tested.
Groups by (sensitivity, specificity) descending, then test_result_lag ascending,
followed by percent_tested (ascending).

This ensures perfect tests (sens=1.0, spec=1.0) come first, with lag=0 before lag=14,
followed by imperfect tests in descending order of accuracy.

# Arguments
- `test_specifications`: Collection of test specifications from a StructVector
- `percent_tested`: Collection of percent_tested values from a StructVector

# Returns
- Vector of indices representing the sorted order

# Example
```julia
# Reorder a StructVector by test specification and percent_tested
results = # ... StructVector{OptimizationResult}
sorted_indices = sort_structvector_by_test_specification_and_percent(
    results.test_specification,
    results.percent_tested
)
sorted_results = results[sorted_indices]
```
"""
function sort_structvector_by_test_specification_and_percent(
        test_specifications,
        percent_tested
    )
    # Sort by (sensitivity, specificity) ascending, then test_result_lag descending,
    # then percent_tested (ascending)
    # This gives: 0.85, 0.9, 1.0 (lag=14), 1.0 (lag=0)
    return sortperm(
        collect(zip(test_specifications, percent_tested));
        by = x -> (
            x[1].sensitivity,  # Ascending order (0.85, 0.9, 1.0)
            x[1].specificity,  # Ascending order (0.85, 0.9, 1.0)
            -x[1].test_result_lag,  # Descending order (14 before 0)
            x[2],  # percent_tested in ascending order
        )
    )
end


"""
    sort_structvector_by_test_specification(test_specifications)

Sort indices of test specifications for reordering a StructVector.
Uses the same ordering as `sort_test_specifications`: test_result_lag (descending),
then (sensitivity, specificity) (ascending).

This function differs from `sort_test_specifications` in that it returns the sorted
**indices** rather than the sorted collection itself, which is needed for reordering
StructVectors while maintaining their structure.

# Arguments
- `test_specifications`: Collection of test specifications from a StructVector

# Returns
- Vector of indices representing the sorted order

# Example
```julia
# Reorder a StructVector by test specification
results = # ... StructVector{OptimizationResult}
sorted_indices = sort_structvector_by_test_specification(results.test_specification)
sorted_results = results[sorted_indices]
```
"""
function sort_structvector_by_test_specification(test_specifications)
    # Use sortperm with the same comparison logic as sort_test_specifications
    # First sort by test_result_lag (descending), then by (sensitivity, specificity) (ascending)

    # First sort those indices by test_result_lag (descending)
    indices_by_lag = sortperm(
        test_specifications;
        by = t -> t.test_result_lag,
        rev = true
    )

    # Then indices sorted by (sensitivity, specificity)
    final_indices = sortperm(
        test_specifications[indices_by_lag];
        by = t -> (t.sensitivity, t.specificity),
        rev = false

    )

    # Map back to original indices
    return indices_by_lag[final_indices]
end

"""
    get_unique_test_specifications_in_sorted_order(
        test_specifications
    )

Extract unique test specifications in the canonical sorted order.
Groups by (sensitivity, specificity) first (descending), then sorts by test_result_lag
(descending) within each group.

This ensures perfect tests (sens=1.0, spec=1.0) come first, with lag=0 before lag=14,
followed by imperfect tests in descending order of accuracy.

# Arguments
- `test_specifications`: Collection of test specifications (may contain duplicates)

# Returns
- Vector of unique test specifications in sorted order

# Example
```julia
results = # ... StructVector{OptimizationResult}
unique_specs = get_unique_test_specifications_in_sorted_order(
    results.test_specification
)
```
"""
function get_unique_test_specifications_in_sorted_order(
        test_specifications
    )
    unique_specs = unique(test_specifications)

    # Sort by (sensitivity, specificity) ascending, then test_result_lag descending
    # This gives: 0.85, 0.9, 1.0 (lag=14), 1.0 (lag=0)
    # Which when reversed in the legend gives: 1.0 (lag=0), 1.0 (lag=14), 0.9, 0.85
    return sort(
        unique_specs;
        by = spec -> (
            spec.sensitivity,  # Ascending order (0.85, 0.9, 1.0)
            spec.specificity,  # Ascending order (0.85, 0.9, 1.0)
            -spec.test_result_lag,  # Descending order (14 before 0)
        )
    )
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
    available_metrics = [
        :accuracies,
        :detection_delays,
        :proportion_timeseries_in_alert,
        :proportion_timeseries_in_outbreak,
        :proportion_alerts_correct,
        :proportion_outbreaks_detected,
        :unavoidable_cases,
        :alert_durations,
        :outbreak_durations,
        :n_alerts,
        :n_outbreaks,
    ]
    @assert outcome in available_metrics "Unknown outcome: :$outcome. Available metrics: $available_metrics"
    return getproperty(result, outcome)
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
