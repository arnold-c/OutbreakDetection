export calculate_simulation_accuracy

"""
    calculate_simulation_accuracy(alert_vec, outbreak_bounds, accuracy_metric)

Calculate accuracy metrics for a single simulation using pre-computed outbreak bounds.

# Arguments
- `alert_vec`: Binary vector indicating alert status
- `outbreak_bounds`: Matrix with outbreak bounds [start, end, duration, size]
- `accuracy_metric`: AccuracyMetric (BalancedAccuracy or F1)

# Returns
- `accuracy`: Overall accuracy score
- `prop_detected`: Proportion of outbreaks detected (sensitivity)
- `prop_correct`: Proportion of alerts that are correct (PPV)
- `delays`: Vector of detection delays
"""
function calculate_simulation_accuracy(
        alert_vec::Union{Vector{Bool}, BitVector},
        outbreak_thresholds::OutbreakThresholds,
        accuracy_metric::AccuracyMetric,
    )
    # Get alert bounds using RLE
    alert_rle = StatsBase.rle(alert_vec)
    alert_bounds = calculate_above_threshold_bounds(alert_rle)

    matched_outbreak_thresholds = match_outbreak_detection_bounds(
        outbreak_thresholds,
        alert_bounds
    )

    # Calculate accuracy based on metric
    accuracy = calculate_accuracy(
        accuracy_metric,
        calculate_ppv(matched_outbreak_thresholds),
        calculate_sensitivity(matched_outbreak_thresholds),
    )

    return accuracy
end
