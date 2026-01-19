export create_test_positive_container

"""
    create_test_positive_container(alert_method, test_positives_vec)

Create a StructVector of TestPositiveContainer objects for ensemble simulations.

This is the main entry point for creating test positive containers from raw test
positive data across multiple simulations. It dispatches to specialized methods
based on the alert method type and returns a StructVector where each element
contains the processed test positive data for a single simulation.

# Arguments
- `alert_method::AlertMethod`: Alert detection method (DailyThreshold, MovingAverage,
  or DailyThresholdMovingAverage)
- `test_positives_vec::Vector{Vector{Int64}}`: Vector of test positive counts for
  each simulation, where each inner vector contains daily counts for one simulation

# Returns
- `StructVector{TestPositiveContainer}`: StructVector containing one TestPositiveContainer
  per simulation, with data processed according to the alert method

# Implementation Notes
This function uses a StructVector approach for memory efficiency and type stability.
Each simulation's data is stored in a separate container, allowing for efficient
indexing and parallel processing during ensemble calculations.

# Examples
```julia
# Create containers for DailyThreshold method
test_positives = [[0, 5, 10, 8, 3], [0, 3, 8, 6, 2]]
alert_method = AlertMethod(DailyThreshold())
containers = create_test_positive_container(alert_method, test_positives)

# Create containers for MovingAverage method
alert_method = AlertMethod(MovingAverage(window = 7))
containers = create_test_positive_container(alert_method, test_positives)
```

# See also
- [`AlertMethod`](@ref)
- [`TestPositiveContainer`](@ref)
- [`SingleAlertTestPositiveContainer`](@ref)
- [`DualAlertTestPositiveContainer`](@ref)
"""
function create_test_positive_container(
        alert_method::AlertMethod,
        test_positives_vec::Vector{Vector{Int64}}
    )
    inner_containers = [
        create_test_positive_container(
                LightSumTypes.variant(alert_method),
                test_positives
            ) for test_positives in test_positives_vec
    ]
    return StructArrays.StructVector(inner_containers)
end

"""
    create_test_positive_container(alert_method::MovingAverage, test_positives)

Create a TestPositiveContainer for MovingAverage alert method.

This method processes test positive data by calculating the moving average over
the specified window. The resulting container stores only the moving average
values, as these are the only values needed for alert generation with this method.

# Arguments
- `alert_method::MovingAverage`: MovingAverage alert method with window parameter
- `test_positives::Vector{Int64}`: Daily test positive counts for a single simulation

# Returns
- `TestPositiveContainer`: Container wrapping a SingleAlertTestPositiveContainer
  with moving average values

# Implementation Notes
The moving average is calculated using a backward-looking window. For days where
the full window is not available (early in the time series), the average is
calculated over the available days.

# See also
- [`MovingAverage`](@ref)
- [`calculate_movingavg`](@ref)
- [`SingleAlertTestPositiveContainer`](@ref)
"""
function create_test_positive_container(
        alert_method::MovingAverage,
        test_positives::Vector{Int64}
    )
    movingavg_test_positives = calculate_movingavg(
        test_positives,
        alert_method.window
    )
    return TestPositiveContainer(
        SingleAlertTestPositiveContainer(
            test_positive_results = movingavg_test_positives
        )
    )
end

"""
    create_test_positive_container(alert_method::DailyThreshold, test_positives)

Create a TestPositiveContainer for DailyThreshold alert method.

This method creates a container that stores the raw daily test positive counts,
as no preprocessing is needed for the DailyThreshold alert method.

# Arguments
- `alert_method::DailyThreshold`: DailyThreshold alert method
- `test_positives::Vector{Int64}`: Daily test positive counts for a single simulation

# Returns
- `TestPositiveContainer`: Container wrapping a SingleAlertTestPositiveContainer
  with raw daily counts

# See also
- [`DailyThreshold`](@ref)
- [`SingleAlertTestPositiveContainer`](@ref)
"""
function create_test_positive_container(
        alert_method::DailyThreshold,
        test_positives::Vector{Int64}
    )
    return TestPositiveContainer(
        SingleAlertTestPositiveContainer(
            test_positive_results = test_positives
        )
    )
end

"""
    create_test_positive_container(alert_method::DailyThresholdMovingAverage, test_positives)

Create a TestPositiveContainer for DailyThresholdMovingAverage alert method.

This method creates a container that stores both the raw daily test positive
counts and their moving averages. Both are needed because this alert method
triggers when either the daily count or the moving average exceeds the threshold.

# Arguments
- `alert_method::DailyThresholdMovingAverage`: DailyThresholdMovingAverage alert
  method with window parameter
- `test_positives::Vector{Int64}`: Daily test positive counts for a single simulation

# Returns
- `TestPositiveContainer`: Container wrapping a DualAlertTestPositiveContainer
  with both raw daily counts and moving average values

# Implementation Notes
The moving average is calculated using a backward-looking window. For days where
the full window is not available (early in the time series), the average is
calculated over the available days.

# See also
- [`DailyThresholdMovingAverage`](@ref)
- [`calculate_movingavg`](@ref)
- [`DualAlertTestPositiveContainer`](@ref)
"""
function create_test_positive_container(
        alert_method::DailyThresholdMovingAverage,
        test_positives::Vector{Int64}
    )
    movingavg_test_positives = calculate_movingavg(
        test_positives,
        alert_method.window
    )
    return TestPositiveContainer(
        DualAlertTestPositiveContainer(
            test_positives = test_positives,
            movingavg_test_positives = movingavg_test_positives
        )
    )
end
