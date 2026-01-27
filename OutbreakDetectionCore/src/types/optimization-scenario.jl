export ScenarioSpecificationVecs,
    OptimizationScenario


"""
    ScenarioSpecificationVecs

Container for vectors of scenario specification parameters used in optimization.

This struct holds vectors of all possible values for each scenario parameter dimension.
These vectors are used to generate the Cartesian product of all parameter combinations,
creating a comprehensive set of optimization scenarios for threshold optimization studies.

# Fields

  - `ensemble_specification_vec::Vector{EnsembleSpecification}`: Disease ensemble
    configurations defining epidemiological dynamics
  - `noise_level_vec::Vector{Float64}`: Observation noise intensity levels to test
  - `noise_type_description_vec::Vector{Symbol}`: Noise types (`:static` for
    time-invariant noise, `:dynamic` for time-varying noise)
  - `test_specification_vec::Vector{IndividualTestSpecification}`: Diagnostic test
    characteristics (sensitivity, specificity)
  - `percent_tested_vec::Vector{Float64}`: Population testing coverage fractions
  - `alert_method_vec::Vector{AlertMethod}`: Alert generation methods (e.g., moving
    average, EWMA)
  - `accuracy_metric_vec::Vector{AccuracyMetric}`: Performance metrics for threshold
    optimization (e.g., F1 score, Matthews correlation)
  - `threshold_bounds_vec::Vector{@NamedTuple{lower::Float64, upper::Float64}}`:
    Search space bounds for threshold optimization
  - `outbreak_specification_vec::Vector{OutbreakSpecification}`: Outbreak detection
    criteria and definitions
  - `alert_filtering_strategy_vec::Vector{AlertFilteringStrategy}`: Strategies for
    filtering alerts during outbreak matching (e.g., all alerts or only post-outbreak)

# Example

```julia
spec_vecs = ScenarioSpecificationVecs(;
    ensemble_specification_vec = [measles_ensemble],
    noise_level_vec = [0.1, 0.2, 0.3],
    noise_type_description_vec = [:static, :dynamic],
    test_specification_vec = [perfect_test, imperfect_test],
    percent_tested_vec = [0.5, 0.8],
    alert_method_vec = [MovingAverage(7)],
    accuracy_metric_vec = [F1Score()],
    threshold_bounds_vec = [(lower = 0.0, upper = 1.0)],
    outbreak_specification_vec = [outbreak_spec],
    alert_filtering_strategy_vec = [AlertFilteringStrategy(AllAlerts())],
)
```

# See Also

  - [`create_scenarios_structvector`](@ref): Function that uses this struct to
    generate scenario combinations
  - [`OptimizationScenario`](@ref): Individual scenario struct created from these
    vectors
  - [`AlertFilteringStrategy`](@ref): Sum type for alert filtering strategies
"""
Base.@kwdef struct ScenarioSpecificationVecs
    ensemble_specification_vec::Vector{EnsembleSpecification}
    noise_level_vec::Vector{Float64}
    noise_type_description_vec::Vector{Symbol}
    test_specification_vec::Vector{IndividualTestSpecification}
    percent_tested_vec::Vector{Float64}
    alert_method_vec::Vector{AlertMethod}
    accuracy_metric_vec::Vector{AccuracyMetric}
    threshold_bounds_vec::Vector{@NamedTuple{lower::Float64, upper::Float64}}
    outbreak_specification_vec::Vector{OutbreakSpecification}
    alert_filtering_strategy_vec::Vector{AlertFilteringStrategy}
end

"""
    OptimizationScenario

Complete specification for a single threshold optimization scenario.

This struct encapsulates all parameters defining a unique optimization scenario,
including disease dynamics, observation noise, diagnostic testing, alert generation,
and outbreak detection criteria. Each scenario represents one combination of parameters
for which optimal detection thresholds will be determined.

The inner constructor validates that:

  - `noise_type_description` is either `:static` or `:dynamic`
  - `threshold_bounds.lower < threshold_bounds.upper`
  - `threshold_bounds.lower >= 0.0`

# Fields

  - `ensemble_specification::EnsembleSpecification`: Disease ensemble configuration
    defining epidemiological dynamics
  - `noise_level::Float64`: Observation noise intensity level
  - `noise_type_description::Symbol`: Noise type (`:static` for time-invariant,
    `:dynamic` for time-varying)
  - `test_specification::IndividualTestSpecification`: Diagnostic test
    characteristics (sensitivity, specificity)
  - `percent_tested::Float64`: Fraction of population tested
  - `alert_method::AlertMethod`: Alert generation method (e.g., moving average)
  - `accuracy_metric::AccuracyMetric`: Performance metric for optimization
    (e.g., F1 score)
  - `threshold_bounds::@NamedTuple{lower::Float64, upper::Float64}`: Search space
    bounds for threshold optimization
  - `outbreak_specification::OutbreakSpecification`: Outbreak detection criteria
  - `alert_filtering_strategy::AlertFilteringStrategy`: Strategy for filtering
    alerts during outbreak matching (default: `AllAlerts`)

# Example

```julia
scenario = OptimizationScenario(;
    ensemble_specification = measles_ensemble,
    noise_level = 0.2,
    noise_type_description = :dynamic,
    test_specification = IndividualTestSpecification(0.95, 0.99),
    percent_tested = 0.8,
    alert_method = MovingAverage(7),
    accuracy_metric = F1Score(),
    threshold_bounds = (lower = 0.0, upper = 1.0),
    outbreak_specification = outbreak_spec,
    alert_filtering_strategy = AlertFilteringStrategy(AllAlerts()),
)
```

# See Also

  - [`ScenarioSpecificationVecs`](@ref): Container for vectors used to generate
    multiple scenarios
  - [`create_scenarios_structvector`](@ref): Function that creates scenario
    combinations
  - [`AlertFilteringStrategy`](@ref): Sum type for alert filtering strategies
"""
Base.@kwdef struct OptimizationScenario
    ensemble_specification::EnsembleSpecification
    noise_level::Float64
    noise_type_description::Symbol
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    alert_method::AlertMethod
    accuracy_metric::AccuracyMetric
    threshold_bounds::@NamedTuple{lower::Float64, upper::Float64}
    outbreak_specification::OutbreakSpecification
    alert_filtering_strategy::AlertFilteringStrategy

    function OptimizationScenario(
            ensemble_specification,
            noise_level,
            noise_type_description,
            test_specification,
            percent_tested,
            alert_method,
            accuracy_metric,
            threshold_bounds,
            outbreak_specification,
            alert_filtering_strategy,
        )
        @assert noise_type_description in [:static, :dynamic] "The noise type must be either :static or :dynamic. Received $noise_type_description"
        @assert threshold_bounds.lower < threshold_bounds.upper
        @assert threshold_bounds.lower >= 0.0
        return new(
            ensemble_specification,
            noise_level,
            noise_type_description,
            test_specification,
            percent_tested,
            alert_method,
            accuracy_metric,
            threshold_bounds,
            outbreak_specification,
            alert_filtering_strategy
        )
    end
end
