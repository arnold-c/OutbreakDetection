export get_most_recent_optimization_filepath, load_optimization_results,
    reshape_results_by_parameter, _extract_parameter

"""
    get_most_recent_optimization_filepath(filename_base, filedir)

Get the most recent optimization file path based on datetime in filename.

Returns Try.Ok(filepath) or Try.Err(message).
"""
function get_most_recent_optimization_filepath(
        filename_base,
        filedir,
    )
    @assert isdir(filedir)
    optimization_files = readdir(filedir)

    if length(optimization_files) == 0
        return Try.Err("No optimization files found.")
    end

    filter_regex = Regex("(.*)$(filename_base)\$")

    filtered_optimization_files = filter(
        f -> contains(f, filter_regex),
        optimization_files,
    )

    if length(filtered_optimization_files) == 0
        return Try.Err("No optimization files found.")
    end

    filtered_optimization_datetimes = Vector{Union{Try.Ok, Try.Err}}(
        undef, length(filtered_optimization_files)
    )

    for (i, f) in pairs(filtered_optimization_files)
        matches = match(filter_regex, f)
        if isnothing(matches)
            filtered_optimization_datetimes[i] = Try.Err(
                "No matches for filename $(f)"
            )
            continue
        end

        filtered_optimization_datetimes[i] = Try.Ok(
            tryparse(
                Dates.DateTime,
                strip(
                    matches[1],
                    '_',
                ),
            ),
        )
    end

    filtered_optimization_datetimes = filter(
        Try.isok, filtered_optimization_datetimes
    )

    if length(filtered_optimization_datetimes) == 0
        return Try.Err("No optimization files found.")
    end

    most_recent_optimization_datetime = sort(
        Try.unwrap.(filtered_optimization_datetimes)
    )[end]

    most_recent_filepath = joinpath(
        filedir,
        string(most_recent_optimization_datetime) *
            "_$(filename_base)",
    )
    return Try.Ok(most_recent_filepath)
end

# ============================================================================
# New result loading and reshaping functions for OptimizationResult
# ============================================================================

"""
    load_optimization_results(results_dir)

Load optimization results from disk.

# Arguments
- `results_dir::String`: Directory containing results

# Returns
- `StructVector{OptimizationResult}`: Loaded results

# Examples
```julia
results = load_optimization_results("results/optimization")

# Efficient filtering with StructVector
high_R0 = filter(r -> r.scenario_params.target_dynamics.R_0 > 15.0, results)
```
"""
function load_optimization_results(results_dir::String)
    results_file = joinpath(results_dir, "optimization_results.jld2")

    if !isfile(results_file)
        error("Results file not found: $results_file")
    end

    return JLD2.load(results_file, "results")
end

"""
    reshape_results_by_parameter(results, parameter_name)

Reshape results by a specific parameter for analysis.

# Arguments
- `results::StructVector{OptimizationResult}`: Results to reshape
- `parameter_name::Symbol`: Parameter to group by

# Returns
- `Dict`: Results grouped by parameter value

# Examples
```julia
results = load_optimization_results("results/optimization")

# Group by R_0
by_R0 = reshape_results_by_parameter(results, :R_0)

# Group by noise level
by_noise = reshape_results_by_parameter(results, :noise_level)

# Analyze each group
for (R0_value, group_results) in by_R0
    mean_acc = Statistics.mean(group_results.accuracy)
    println("R_0 = \$R0_value: mean accuracy = \$mean_acc")
end
```
"""
function reshape_results_by_parameter(
        results::StructVector{OptimizationResult}, parameter_name::Symbol
    )
    # Extract parameter values
    param_values = _extract_parameter(results, parameter_name)

    # Group by unique values
    unique_values = unique(param_values)

    grouped = Dict()
    for value in unique_values
        mask = param_values .== value
        grouped[value] = results[mask]
    end

    return grouped
end

"""
    _extract_parameter(results, parameter_name)

Internal function to extract parameter values from results.

Handles nested parameter structures in ScenarioParameters.

# Arguments
- `results::StructVector{OptimizationResult}`: Results to extract from
- `parameter_name::Symbol`: Parameter to extract

# Returns
- `Vector`: Extracted parameter values
"""
function _extract_parameter(
        results::StructVector{OptimizationResult}, parameter_name::Symbol
    )
    # Map common parameter names to their locations in ScenarioParameters
    if parameter_name == :R_0
        return [r.scenario_params.target_dynamics.R_0 for r in results]
    elseif parameter_name == :latent_period_days
        return [r.scenario_params.target_dynamics.latent_period_days for r in results]
    elseif parameter_name == :infectious_duration_days
        return [
            r.scenario_params.target_dynamics.infectious_duration_days for r in results
        ]
    elseif parameter_name == :test_sensitivity
        return [r.scenario_params.test_spec.sensitivity for r in results]
    elseif parameter_name == :test_specificity
        return [r.scenario_params.test_spec.specificity for r in results]
    elseif parameter_name == :percent_clinic_tested
        return [r.scenario_params.percent_clinic_tested for r in results]
    elseif parameter_name == :outbreak_threshold
        return [r.scenario_params.outbreak_threshold for r in results]
    else
        error("Unknown parameter name: $parameter_name")
    end
end
