export calculate_mean_ensemble_accuracy

"""
    calculate_mean_ensemble_accuracy(
        alert_threshold,
        accuracy_metric,
        test_positives_container,
        outbreak_status_vecs,
        outbreak_bounds_vecs
    )

Calculate accuracy metrics across all simulations in the ensemble.

# Arguments
- `alert_threshold`: Threshold for generating alerts
- `accuracy_metric`: AccuracyMetric for evaluation
- `test_positives_container`: TestPositiveContainer with test data
- `outbreak_status_vecs`: Vector of outbreak status vectors
- `outbreak_bounds_vecs`: Vector of outbreak bounds matrices

# Returns
- `mean_accuracy`: Mean accuracy across simulations
"""
function calculate_mean_ensemble_accuracy(
        alert_threshold::Float64,
        accuracy_metric::AccuracyMetric,
        test_positives_container::StructVector{TestPositiveContainer},
        outbreak_thresholds::StructVector{OutbreakThresholds},
    )
    nsims = length(outbreak_thresholds)

    # Accumulators for metrics
    accuracies = Vector{Float64}(undef, nsims)

    for sim in eachindex(accuracies)
        # Generate alerts based on threshold and alert method
        alert_vec = generate_alerts(
            test_positives_container[sim],
            alert_threshold
        )

        # Get pre-computed outbreak bounds for this simulation
        outbreak_bounds = outbreak_thresholds[sim]

        # Calculate detection characteristics
        sim_accuracy = calculate_simulation_accuracy(
            alert_vec,
            outbreak_bounds,
            accuracy_metric
        )

        accuracies[sim] = sim_accuracy
    end

    mean_accuracy = NaNMath.mean(accuracies)

    return mean_accuracy
end
