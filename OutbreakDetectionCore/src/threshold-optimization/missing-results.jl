export find_missing_scenarios

"""
    find_missing_scenarios(
		all_scenarios::StructVector{OptimizationScenario},
        completed_results::StructVector{OptimizationResult}
	) -> StructVector{OptimizationScenario}

Find optimization scenarios that haven't been computed yet by comparing against completed results.

This function efficiently identifies which scenarios from a complete set still need to be
optimized by comparing them against previously completed results. It uses a dictionary-based
lookup for O(1) comparison performance, making it suitable for large scenario sets.

# Arguments
- `all_scenarios::StructVector{OptimizationScenario}`: Complete set of scenarios to check
- `completed_results::StructVector{OptimizationResult}`: Previously completed optimization results

# Returns
- `StructVector{OptimizationScenario}`: Subset of scenarios that haven't been computed yet.
  Returns all scenarios if `completed_results` is empty.

# Algorithm
1. If no completed results exist, returns all scenarios
2. Builds a dictionary of completed scenario keys for O(1) lookup
3. Filters scenarios by checking if their key exists in the completed dictionary
4. Returns only scenarios not found in completed results

# Example
```julia
# Create scenarios and run some optimizations
all_scenarios = create_scenarios_structvector(specification_vecs)
completed_results = load_previous_optimization_results_structvector(results_dir)

# Find which scenarios still need to be computed
missing = find_missing_scenarios(all_scenarios, completed_results)
# Returns only scenarios not in completed_results

# Run optimizations on missing scenarios
for scenario in missing
    result = threshold_optimization(scenario, ...)
    save_result(result)
end
```

# Performance
- Time complexity: O(n + m) where n = number of scenarios, m = number of results
- Space complexity: O(m) for the completed scenarios dictionary

# See Also
- [`scenario_key`](@ref): Generates keys for scenario comparison
- [`result_key`](@ref): Generates keys from optimization results
- [`build_scenario_dict`](@ref): Creates the lookup dictionary
"""
function find_missing_scenarios(
        all_scenarios::StructVector{OptimizationScenario},
        completed_results::StructVector{OptimizationResult}
    )
    if isempty(completed_results)
        return all_scenarios
    end

    # Build a dictionary of completed scenarios for O(1) lookup
    completed_dict = build_scenario_dict(completed_results)

    # Filter to find missing scenarios using dictionary lookup
    missing_mask = map(all_scenarios) do scenario
        !haskey(completed_dict, scenario_key(scenario))
    end

    return all_scenarios[missing_mask]
end

"""
    build_scenario_dict(results::StructVector{OptimizationResult}) -> Dict{Tuple, Bool}

Build a dictionary mapping scenario keys to true for all completed results.

This function creates an efficient lookup structure for checking whether a scenario
has already been computed. The dictionary uses scenario keys (tuples of scenario
parameters) as keys and `true` as values for O(1) lookup performance.

# Arguments
- `results::StructVector{OptimizationResult}`: Completed optimization results

# Returns
- `Dict{Tuple, Bool}`: Dictionary mapping scenario keys to `true` for all completed results

# See Also
- [`result_key`](@ref): Generates keys from optimization results
- [`scenario_key`](@ref): Generates keys from optimization scenarios
- [`find_missing_scenarios`](@ref): Uses this dictionary to identify missing scenarios
"""
function build_scenario_dict(
        results::StructVector{OptimizationResult},
    )::Dict{Tuple, Bool}
    dict = Dict{Tuple, Bool}()

    for result in results
        dict[result_key(result)] = true
    end
    return dict
end

"""
    scenario_key(scenario::OptimizationScenario) -> Tuple

Generate a unique tuple key from an `OptimizationScenario` for identification and comparison.

This function extracts all relevant parameters from an optimization scenario and
combines them into a tuple that uniquely identifies the scenario. The key is used
for efficient lookup and comparison operations when identifying missing scenarios.

# Arguments
- `scenario::OptimizationScenario`: The optimization scenario to generate a key for

# Returns
- `Tuple`: A tuple containing all scenario parameters in the following order:
  - `ensemble_specification`: Ensemble configuration
  - `noise_level`: Level of noise in the simulation
  - `noise_type_description`: Description of noise type
  - `test_specification`: Diagnostic testing parameters
  - `percent_tested`: Percentage of population tested
  - `alert_method`: Method used for generating alerts
  - `accuracy_metric`: Metric used to evaluate accuracy
  - `threshold_bounds`: Bounds for threshold optimization
  - `outbreak_specification`: Outbreak parameters

# See Also
- [`result_key`](@ref): Generates equivalent keys from optimization results
- [`build_scenario_dict`](@ref): Uses these keys to build lookup dictionaries
- [`find_missing_scenarios`](@ref): Uses these keys to identify missing scenarios
"""
function scenario_key(scenario::OptimizationScenario)
    return (
        scenario.ensemble_specification,
        scenario.noise_level,
        scenario.noise_type_description,
        scenario.test_specification,
        scenario.percent_tested,
        scenario.alert_method,
        scenario.accuracy_metric,
        scenario.threshold_bounds,
        scenario.outbreak_specification,
    )
end

"""
    result_key(result::OptimizationResult) -> Tuple

Generate a unique tuple key from an `OptimizationResult` for identification and comparison.

This function extracts all relevant parameters from an optimization result and
combines them into a tuple that uniquely identifies the scenario that was optimized.
The key structure matches that of [`scenario_key`](@ref), enabling direct comparison
between completed results and pending scenarios.

# Arguments
- `result::OptimizationResult`: The optimization result to generate a key for

# Returns
- `Tuple`: A tuple containing all scenario parameters in the following order:
  - `ensemble_specification`: Ensemble configuration
  - `noise_level`: Level of noise in the simulation
  - `noise_type_description`: Description of noise type
  - `test_specification`: Diagnostic testing parameters
  - `percent_tested`: Percentage of population tested
  - `alert_method`: Method used for generating alerts
  - `accuracy_metric`: Metric used to evaluate accuracy
  - `threshold_bounds`: Bounds for threshold optimization
  - `outbreak_specification`: Outbreak parameters

# See Also
- [`scenario_key`](@ref): Generates equivalent keys from optimization scenarios
- [`build_scenario_dict`](@ref): Uses these keys to build lookup dictionaries
- [`find_missing_scenarios`](@ref): Uses these keys to identify missing scenarios
"""
function result_key(result::OptimizationResult)
    return (
        result.ensemble_specification,
        result.noise_level,
        result.noise_type_description,
        result.test_specification,
        result.percent_tested,
        result.alert_method,
        result.accuracy_metric,
        result.threshold_bounds,
        result.outbreak_specification,
    )
end
