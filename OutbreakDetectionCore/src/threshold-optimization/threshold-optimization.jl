export threshold_optimization

"""
    threshold_optimization(
        scenario,
        test_positives_container,
        ensemble_simulation,
        outbreak_thresholds,
        opt_params
    ) -> OptimizationResult

Optimize alert threshold for outbreak detection using multistart optimization.

This function finds the optimal alert threshold that maximizes detection accuracy
for a given outbreak detection scenario. It uses multistart optimization with
Sobol sequences to explore the parameter space and NLopt for local optimization.
The optimization minimizes the loss function (1 - accuracy) across all simulations
in the ensemble.

After finding the optimal alert threshold, the scenario is re-run (outside of the
optimization hot-loop) with the optimal threshold to produce a set of summary
metrics that can be compared between scenarios.

# Arguments

  - `scenario::OptimizationScenario`: Complete specification of the optimization problem including ensemble parameters, noise characteristics, test specifications, alert method, accuracy metric, and threshold bounds
  - `test_positives_container::StructVector{TestPositiveContainer}`: Container with test positive data for each simulation, structured according to the alert method
  - `ensemble_simulation::StructVector{SEIRRun}`: Collection of SEIR simulation runs containing incidence and infection trajectories
  - `outbreak_thresholds::StructVector{OutbreakThresholds}`: Pre-computed outbreak threshold bounds for each simulation based on the outbreak specification
  - `opt_params::ThresholdOptimizationParameters`: Optimization algorithm parameters including convergence tolerances (xtol_rel, xtol_abs), maximum evaluations (maxeval), and number of Sobol points for multistart

# Returns

  - `OptimizationResult`: Comprehensive results containing the optimal threshold, accuracy metrics (PPV, sensitivity, arithmetic mean/F1 accuracy score), detection performance (delays, unavoidable cases), and temporal characteristics (alert/outbreak durations and proportions)

# Details

The optimization process:

  1. Creates an objective function that evaluates candidate thresholds by calculating mean ensemble accuracy
  2. Configures multistart optimization using TikTak method with Sobol sequences for global exploration
  3. Uses NLopt local optimization (default: LN_NELDERMEAD) for refinement within each multistart iteration
  4. Tracks the best solution found across all optimization iterations
  5. Calculates comprehensive performance metrics at the optimal threshold

# Example

```julia
# Define optimization scenario
scenario = OptimizationScenario(;
    ensemble_specification = ensemble_spec,
    noise_level = 0.1,
    noise_type_description = :static,
    test_specification = :daily,
    percent_tested = 0.5,
    alert_method = DailyThreshold(),
    accuracy_metric = F1Score(),
    threshold_bounds = (lower = 0.0, upper = 100.0),
    outbreak_specification = outbreak_spec
)

# Run optimization
result = threshold_optimization(
    scenario,
    test_positives_container,
    ensemble_simulation,
    outbreak_thresholds,
    ThresholdOptimizationParameters()
)

println("Optimal threshold: ", result.optimal_threshold)
println("Mean accuracy: ", mean(result.accuracies))
```
"""
function threshold_optimization(
        scenario::OptimizationScenario,
        test_positives_container::StructVector{TestPositiveContainer},
        ensemble_simulation::StructVector{SEIRRun},
        outbreak_thresholds::StructVector{OutbreakThresholds},
        opt_params::ThresholdOptimizationParameters,
    )
    # Create tracker instance for this scenario
    tracker = OptimizationTracker()

    # Create objective function closure that updates tracker
    objective = params -> multistart_objective_function(
        params,
        scenario,
        test_positives_container,
        outbreak_thresholds,
        tracker,
    )

    # Setup multistart problem
    problem = MultistartOptimization.MinimizationProblem(
        objective,
        [scenario.threshold_bounds.lower],
        [scenario.threshold_bounds.upper]
    )

    # Configure local optimization method
    local_method = MultistartOptimization.NLoptLocalMethod(
        opt_params.local_algorithm;
        xtol_rel = opt_params.xtol_rel,
        xtol_abs = opt_params.xtol_abs,
        maxeval = opt_params.maxeval,
    )

    # Configure multistart method (TikTak uses Sobol sequences)
    multistart_method = MultistartOptimization.TikTak(opt_params.n_sobol_points)

    # Run optimization
    MultistartOptimization.multistart_minimization(
        multistart_method,
        local_method,
        problem
    )

    optimal_results = calculate_optimal_results(
        tracker.optimal_threshold,
        scenario.accuracy_metric,
        test_positives_container,
        ensemble_simulation,
        outbreak_thresholds,
    )

    return OptimizationResult(
        # From scenario
        ensemble_specification = scenario.ensemble_specification,
        noise_level = scenario.noise_level,
        noise_type_description = scenario.noise_type_description,
        test_specification = scenario.test_specification,
        percent_tested = scenario.percent_tested,
        alert_method = scenario.alert_method,
        accuracy_metric = scenario.accuracy_metric,
        threshold_bounds = scenario.threshold_bounds,
        outbreak_specification = scenario.outbreak_specification,
        # From optimized values
        n_alerts = optimal_results.n_alerts,
        n_outbreaks = optimal_results.n_outbreaks,
        optimal_threshold = tracker.optimal_threshold,
        accuracies = optimal_results.accuracies,
        proportion_alerts_correct = optimal_results.proportion_alerts_correct,
        proportion_outbreaks_detected = optimal_results.proportion_outbreaks_detected,
        detection_delays = optimal_results.detection_delays,
        unavoidable_cases = optimal_results.unavoidable_cases,
        alert_durations = optimal_results.alert_durations,
        outbreak_durations = optimal_results.outbreak_durations,
        proportion_timeseries_in_alert = optimal_results.proportion_timeseries_in_alert,
        proportion_timeseries_in_outbreak = optimal_results.proportion_timeseries_in_outbreak,
    )
end
