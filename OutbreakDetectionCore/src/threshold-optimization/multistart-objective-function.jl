export multistart_objective_function

"""
    multistart_objective_function(
        params,
        scenario,
        test_positives_container,
        outbreak_status_vecs,
        outbreak_bounds_vecs,
        tracker
    )

Objective function for threshold optimization using multistart optimization.

Evaluates a candidate alert threshold by calculating detection accuracy across
all simulations in the ensemble. Updates the tracker with the best solution found.

# Arguments
- `params`: Vector containing the alert threshold parameter
- `scenario`: OptimizationScenario with ensemble and detection specifications
- `test_positives_container`: TestPositiveContainer with test positive data
- `outbreak_status_vecs`: Vector of outbreak status vectors for each simulation
- `outbreak_bounds_vecs`: Vector of outbreak bounds matrices for each simulation
- `tracker`: OptimizationTracker to track best solution

# Returns
- `loss`: Loss value (1 - accuracy) for minimization
"""
function multistart_objective_function(
        params::Vector{Float64},
        scenario::OptimizationScenario,
        test_positives_container::StructVector{TestPositiveContainer},
        outbreak_thresholds::StructVector{OutbreakThresholds},
        tracker::OptimizationTracker,
    )
    # Extract threshold from params (single parameter optimization)
    alert_threshold = params[1]

    # Calculate accuracy across all simulations
    accuracy = calculate_mean_ensemble_accuracy(
        alert_threshold,
        scenario.accuracy_metric,
        test_positives_container,
        outbreak_thresholds,
    )

    # Calculate loss (for minimization)
    loss = 1.0 - accuracy

    # Update tracker if this is the best solution found
    if loss < tracker.best_loss
        tracker.best_loss = loss
        tracker.optimal_threshold = alert_threshold
        tracker.best_accuracy = accuracy
    end

    return loss
end
