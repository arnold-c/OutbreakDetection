export OptimizationResult, ThresholdOptimizationResult

using AutoHashEquals: AutoHashEquals

"""
    OptimizationResult

Base result from an optimization run.

Designed for StructVector storage with AutoHashEquals for efficient lookup.

# Fields
- `scenario_params::ScenarioParameters`: Scenario that was optimized
- `optimal_threshold::Float64`: Optimal detection threshold
- `accuracy::Float64`: Detection accuracy at optimal threshold
- `sensitivity::Float64`: True positive rate
- `specificity::Float64`: True negative rate
- `mean_detection_delay::Float64`: Mean delay to detection (days)
- `proportion_detected::Float64`: Proportion of outbreaks detected

# Examples
```julia
result = OptimizationResult(
    scenario_params = scenario,
    optimal_threshold = 0.05,
    accuracy = 0.95,
    sensitivity = 0.92,
    specificity = 0.98,
    mean_detection_delay = 14.5,
    proportion_detected = 0.88
)

# Store in StructVector
results = StructVector{OptimizationResult}([result1, result2, ...])

# Efficient filtering
high_accuracy = filter(r -> r.accuracy > 0.9, results)
```
"""
AutoHashEquals.@auto_hash_equals struct OptimizationResult
    scenario_params::ScenarioParameters
    optimal_threshold::Float64
    accuracy::Float64
    sensitivity::Float64
    specificity::Float64
    mean_detection_delay::Float64
    proportion_detected::Float64
end

"""
    ThresholdOptimizationResult

Extended result with detailed threshold characteristics.

Includes full threshold sweep results for analysis.

# Fields
- `base_result::OptimizationResult`: Base optimization result
- `threshold_sweep::Vector{Float64}`: Thresholds tested
- `accuracy_sweep::Vector{Float64}`: Accuracy at each threshold
- `sensitivity_sweep::Vector{Float64}`: Sensitivity at each threshold
- `specificity_sweep::Vector{Float64}`: Specificity at each threshold

# Examples
```julia
result = ThresholdOptimizationResult(
    base_result = opt_result,
    threshold_sweep = [0.01, 0.02, 0.05, 0.1],
    accuracy_sweep = [0.85, 0.92, 0.95, 0.88],
    sensitivity_sweep = [0.95, 0.92, 0.88, 0.75],
    specificity_sweep = [0.75, 0.92, 0.98, 0.99]
)

# Access optimal values
result.base_result.optimal_threshold  # 0.05
result.base_result.accuracy  # 0.95

# Analyze threshold sweep
using Plots
plot(result.threshold_sweep, result.accuracy_sweep, label="Accuracy")
```
"""
AutoHashEquals.@auto_hash_equals struct ThresholdOptimizationResult
    base_result::OptimizationResult
    threshold_sweep::Vector{Float64}
    accuracy_sweep::Vector{Float64}
    sensitivity_sweep::Vector{Float64}
    specificity_sweep::Vector{Float64}
end
