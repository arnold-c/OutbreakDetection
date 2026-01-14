export create_test_positive_container

function create_test_positive_container(
        alert_method::AlertMethod,
        test_positives_vec::Vector{Vector{Int64}}
    )
    return TestPositiveContainer(
        create_test_positive_container(
            LightSumTypes.variant(alert_method),
            test_positives_vec
        )
    )
end

function create_test_positive_container(
        alert_method::MovingAverage,
        test_positives_vec
    )
    movingavg_test_positives = calculate_movingavg(
        test_positives_vec,
        alert_method.window
    )
    return SingleAlertTestPositiveContainer(
        test_positive_results = movingavg_test_positives
    )
end

function create_test_positive_container(
        alert_method::DailyThreshold,
        test_positives_vec
    )
    return SingleAlertTestPositiveContainer(
        test_positive_results = test_positives_vec
    )
end

function create_test_positive_container(
        alert_method::DailyThresholdMovingAverage,
        test_positives_vec
    )
    movingavg_test_positives = calculate_movingavg(
        test_positives_vec,
        alert_method.window
    )
    return DualAlertTestPositiveContainer(
        test_positives = test_positives_vec,
        movingavg_test_positives = movingavg_test_positives
    )
end
