export generate_alerts

"""
    generate_alerts(container, threshold)

Generate binary alert vector based on threshold and alert method.

# Arguments
- `container`: TestPositiveContainer with test positive data for a single simulation
- `threshold`: Alert threshold value

# Returns
- Binary vector indicating alert status at each time point
"""
function generate_alerts(
        container::TestPositiveContainer,
        threshold::Float64,
    )
    return _generate_alerts(
        LightSumTypes.variant(container),
        threshold
    )
end

# Dispatch for DailyThreshold or MovingAverage (SingleAlertTestPositiveContainers)
function _generate_alerts(
        container::SingleAlertTestPositiveContainer,
        threshold::Float64,
    )
    test_positives = container.test_positive_results
    return vec(test_positives .>= threshold)
end

# Dispatch for DailyThresholdMovingAverage (DualAlertTestPositiveContainer)
function _generate_alerts(
        container::DualAlertTestPositiveContainer,
        threshold::Float64,
    )
    test_positives = container.test_positives
    movingavg = container.movingavg_test_positives

    # Alert if EITHER exceeds threshold
    return vec((test_positives .>= threshold) .| (movingavg .>= threshold))
end
