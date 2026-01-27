export create_scenarios_structvector

"""
    create_scenarios_structvector(
    	specification_vecs::ScenarioSpecificationVecs
	) -> StructVector{OptimizationScenario}

Create a StructVector of optimization scenarios from specification vectors.

This function generates all possible combinations of scenario parameters by taking
the Cartesian product of the input specification vectors. Each combination represents
a unique optimization scenario with specific ensemble, noise, testing, alert, and
outbreak characteristics.

The function performs the following steps:

  1. Unpacks all specification vectors from the input struct
  2. Computes the Cartesian product of all specification vectors
  3. Creates an `OptimizationScenario` for each combination
  4. Returns a `StructVector` for efficient columnar storage and access

# Arguments

  - `specification_vecs::ScenarioSpecificationVecs`: A struct containing vectors of
    specifications for all scenario parameters:
      - `ensemble_specification_vec`: Disease ensemble configurations
      - `noise_level_vec`: Noise intensity levels
      - `noise_type_description_vec`: Noise types (`:static` or `:dynamic`)
      - `test_specification_vec`: Diagnostic test specifications
      - `percent_tested_vec`: Testing coverage percentages
      - `alert_method_vec`: Alert generation methods
      - `accuracy_metric_vec`: Performance metrics for optimization
      - `threshold_bounds_vec`: Search bounds for threshold optimization
      - `outbreak_specification_vec`: Outbreak detection criteria
      - `alert_filtering_strategy_vec`: Alert filtering strategies

# Returns

  - `StructVector{OptimizationScenario}`: A columnar storage structure containing all
    scenario combinations, enabling efficient vectorized operations and memory access

# Example

```julia
# Define specification vectors
spec_vecs = ScenarioSpecificationVecs(;
    ensemble_specification_vec = [measles_ensemble, flu_ensemble],
    noise_level_vec = [0.1, 0.2],
    noise_type_description_vec = [:static, :dynamic],
    test_specification_vec = [test_spec_1],
    percent_tested_vec = [0.5, 0.8],
    alert_method_vec = [MovingAverage(7)],
    accuracy_metric_vec = [F1Score()],
    threshold_bounds_vec = [(lower = 0.0, upper = 1.0)],
    outbreak_specification_vec = [outbreak_spec],
    alert_filtering_strategy_vec = [AlertFilteringStrategy(AllAlerts())],
)

# Create all scenario combinations
scenarios = create_scenarios_structvector(spec_vecs)
# Returns 2 × 2 × 2 × 1 × 2 × 1 × 1 × 1 × 1 × 1 = 16 scenarios
```

# See Also

  - [`ScenarioSpecificationVecs`](@ref): Input struct containing specification vectors
  - [`OptimizationScenario`](@ref): Individual scenario struct
  - [`StructArrays.StructVector`](@ref): Columnar storage format for scenarios
"""
function create_scenarios_structvector(specification_vecs::ScenarioSpecificationVecs)
    UnPack.@unpack ensemble_specification_vec,
        noise_level_vec,
        noise_type_description_vec,
        test_specification_vec,
        percent_tested_vec,
        alert_method_vec,
        accuracy_metric_vec,
        threshold_bounds_vec,
        outbreak_specification_vec,
        alert_filtering_strategy_vec = specification_vecs

    combinations = Iterators.product(
        ensemble_specification_vec,
        noise_level_vec,
        noise_type_description_vec,
        test_specification_vec,
        percent_tested_vec,
        alert_method_vec,
        accuracy_metric_vec,
        threshold_bounds_vec,
        outbreak_specification_vec,
        alert_filtering_strategy_vec,
    )
    n_combinations = length(combinations)

    scenarios_vec = Vector{OptimizationScenario}(undef, n_combinations)

    for (
            i, (
                ensemble_spec,
                noise_level,
                noise_type_description,
                test_spec,
                percent_tested,
                alert_method,
                accuracy_metric,
                threshold_bounds,
                outbreak_spec,
                alert_filtering_strategy,
            ),
        ) in enumerate(combinations)

        scenarios_vec[i] = OptimizationScenario(;
            ensemble_specification = ensemble_spec,
            noise_level = noise_level,
            noise_type_description = noise_type_description,
            test_specification = test_spec,
            percent_tested = percent_tested,
            alert_method = alert_method,
            accuracy_metric = accuracy_metric,
            threshold_bounds = threshold_bounds,
            outbreak_specification = outbreak_spec,
            alert_filtering_strategy = alert_filtering_strategy,
        )
    end

    return StructVector(scenarios_vec)
end
