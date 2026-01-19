export TestPositiveContainer,
    SingleAlertTestPositiveContainer,
    DualAlertTestPositiveContainer

"""
    AbstractTestPositiveContainer

Abstract base type for all test positive container variants.

This serves as the supertype for the `TestPositiveContainer` sum type and its
variants. Test positive containers store processed diagnostic test results in
formats optimized for different alert detection methods.

# See also
- [`TestPositiveContainer`](@ref)
- [`SingleAlertTestPositiveContainer`](@ref)
- [`DualAlertTestPositiveContainer`](@ref)
"""
abstract type AbstractTestPositiveContainer end

"""
    SingleAlertTestPositiveContainer{T <: Union{Int64, Float64}} <: AbstractTestPositiveContainer

Container for test positive data used with single-criterion alert methods.

This container stores a single vector of test positive results, used by alert
methods that evaluate only one criterion for triggering alerts. It is used by
both `DailyThreshold` (which stores raw daily counts) and `MovingAverage`
(which stores moving average values).

# Type Parameters
- `T <: Union{Int64, Float64}`: Element type of the test positive results
  - `Int64` for raw daily counts (DailyThreshold method)
  - `Float64` for moving average values (MovingAverage method)

# Fields
- `test_positive_results::Vector{T}`: Time series of test positive values,
  either raw daily counts or moving averages depending on the alert method

# Constructor
    SingleAlertTestPositiveContainer(; test_positive_results)

# Examples
```julia
# Container for DailyThreshold method (raw counts)
daily_container = SingleAlertTestPositiveContainer(
    test_positive_results = [0, 5, 10, 8, 3]
)

# Container for MovingAverage method (moving averages)
movavg_container = SingleAlertTestPositiveContainer(
    test_positive_results = [0.0, 2.5, 5.0, 7.67, 6.5]
)
```

# See also
- [`AbstractTestPositiveContainer`](@ref): Parent abstract type
- [`TestPositiveContainer`](@ref): Sum type wrapper
- [`DualAlertTestPositiveContainer`](@ref): Container for dual-criterion methods
- [`DailyThreshold`](@ref): Alert method using raw daily counts
- [`MovingAverage`](@ref): Alert method using moving averages
"""
Base.@kwdef struct SingleAlertTestPositiveContainer{T <: Union{Int64, Float64}} <: AbstractTestPositiveContainer
    test_positive_results::Vector{T}
end

"""
    DualAlertTestPositiveContainer <: AbstractTestPositiveContainer

Container for test positive data used with dual-criterion alert methods.

This container stores both raw daily test positive counts and their moving
averages, used by the `DailyThresholdMovingAverage` alert method which triggers
when either the daily count or the moving average exceeds the threshold.

# Fields
- `test_positives::Vector{Int64}`: Time series of raw daily test positive counts
- `movingavg_test_positives::Vector{Float64}`: Time series of moving average
  test positive values calculated over a specified window

# Constructor
    DualAlertTestPositiveContainer(; test_positives, movingavg_test_positives)

# Examples
```julia
# Container for DailyThresholdMovingAverage method
dual_container = DualAlertTestPositiveContainer(
    test_positives = [0, 5, 10, 8, 3],
    movingavg_test_positives = [0.0, 2.5, 5.0, 7.67, 6.5]
)

# Access both data streams
daily_values = dual_container.test_positives
movavg_values = dual_container.movingavg_test_positives
```

# Implementation Notes
Both vectors should have the same length, representing the same time period.
The moving average is typically calculated using a backward-looking window,
with partial windows used for early time points.

# See also
- [`AbstractTestPositiveContainer`](@ref): Parent abstract type
- [`TestPositiveContainer`](@ref): Sum type wrapper
- [`SingleAlertTestPositiveContainer`](@ref): Container for single-criterion methods
- [`DailyThresholdMovingAverage`](@ref): Alert method using both criteria
- [`calculate_movingavg`](@ref): Function for computing moving averages
"""
Base.@kwdef struct DualAlertTestPositiveContainer <: AbstractTestPositiveContainer
    test_positives::Vector{Int64}
    movingavg_test_positives::Vector{Float64}
end

"""
    TestPositiveContainer

Sum type for test positive data containers.

This is a LightSumTypes-based sum type that can hold one of two variants:
- `SingleAlertTestPositiveContainer`: For single-criterion alert methods
  (DailyThreshold or MovingAverage)
- `DualAlertTestPositiveContainer`: For dual-criterion alert methods
  (DailyThresholdMovingAverage)

The sum type provides type-safe storage of test positive data in formats
optimized for different alert detection methods, enabling efficient dispatch
and memory usage.

# Variants
- `SingleAlertTestPositiveContainer{Int64}`: Raw daily counts for DailyThreshold
- `SingleAlertTestPositiveContainer{Float64}`: Moving averages for MovingAverage
- `DualAlertTestPositiveContainer`: Both raw counts and moving averages for
  DailyThresholdMovingAverage

# Construction
    TestPositiveContainer(SingleAlertTestPositiveContainer(...))
    TestPositiveContainer(DualAlertTestPositiveContainer(...))

# Examples
```julia
# Create container for DailyThreshold method
daily_container = TestPositiveContainer(
    SingleAlertTestPositiveContainer(test_positive_results = [0, 5, 10, 8, 3])
)

# Create container for MovingAverage method
movavg_container = TestPositiveContainer(
    SingleAlertTestPositiveContainer(test_positive_results = [0.0, 2.5, 5.0])
)

# Create container for DailyThresholdMovingAverage method
dual_container = TestPositiveContainer(
    DualAlertTestPositiveContainer(
        test_positives = [0, 5, 10, 8, 3],
        movingavg_test_positives = [0.0, 2.5, 5.0, 7.67, 6.5]
    )
)

# Pattern matching with LightSumTypes
using LightSumTypes
variant_type = LightSumTypes.variant(daily_container)
```

# Usage in Ensemble Simulations
Typically used within a `StructVector{TestPositiveContainer}` where each element
contains the test positive data for a single simulation in an ensemble.

# See also
- [`AbstractTestPositiveContainer`](@ref): Parent abstract type
- [`SingleAlertTestPositiveContainer`](@ref): Single-criterion variant
- [`DualAlertTestPositiveContainer`](@ref): Dual-criterion variant
- [`create_test_positive_container`](@ref): Factory function for creating containers
- [`AlertMethod`](@ref): Alert method specification
"""
LightSumTypes.@sumtype TestPositiveContainer(SingleAlertTestPositiveContainer, DualAlertTestPositiveContainer) <: AbstractTestPositiveContainer
