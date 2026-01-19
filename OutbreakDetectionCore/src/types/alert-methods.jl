export AlertMethod,
    DailyThreshold,
    MovingAverage,
    DailyThresholdMovingAverage,
    alert_method_to_string

"""
    AbstractAlertMethod

Abstract base type for all alert detection method variants.

This serves as the supertype for the `AlertMethod` sum type and its variants.

# See also
- [`AlertMethod`](@ref)
- [`DailyThreshold`](@ref)
- [`MovingAverage`](@ref)
- [`DailyThresholdMovingAverage`](@ref)
"""
abstract type AbstractAlertMethod end

"""
    DailyThreshold

Alert method variant that triggers alerts based only on daily threshold values.

This method compares daily incidence or test positivity values against a
threshold. An alert is triggered when the daily value meets or exceeds the
threshold.

# See also
- [`AlertMethod`](@ref)
- [`MovingAverage`](@ref)
- [`DailyThresholdMovingAverage`](@ref)
"""
struct DailyThreshold end

"""
    MovingAverage

Alert method variant that triggers alerts based only on moving average values.

This method compares the moving average of incidence or test positivity values
against a threshold. An alert is triggered when the moving average meets or
exceeds the threshold.

# See also
- [`AlertMethod`](@ref)
- [`DailyThreshold`](@ref)
- [`DailyThresholdMovingAverage`](@ref)
"""
Base.@kwdef struct MovingAverage
    window::Int64 = 7
end

"""
    DailyThresholdMovingAverage

Alert method variant that triggers alerts based on both daily and moving
average thresholds.

This method compares both the daily values and moving average values against
a threshold. An alert is triggered when either the daily value or the moving
average meets or exceeds the threshold.

# See also
- [`AlertMethod`](@ref)
- [`DailyThreshold`](@ref)
- [`MovingAverage`](@ref)
"""
Base.@kwdef struct DailyThresholdMovingAverage
    window::Int64 = 7
end

"""
    AlertMethod

Sum type for outbreak alert detection methods.

This is a LightSumTypes-based sum type that can hold one of three variants:
- `DailyThreshold`: Alert based on daily values only
- `MovingAverage`: Alert based on moving average only
- `DailyThresholdMovingAverage`: Alert based on either daily or moving average

# Construction
    AlertMethod(DailyThreshold())
    AlertMethod(MovingAverage())
    AlertMethod(DailyThresholdMovingAverage())

# See also
- [`DailyThreshold`](@ref)
- [`MovingAverage`](@ref)
- [`DailyThresholdMovingAverage`](@ref)
- [`getstring`](@ref)
"""
LightSumTypes.@sumtype AlertMethod(DailyThreshold, MovingAverage, DailyThresholdMovingAverage) <: AbstractAlertMethod

"""
    getstring(alert_method)

Extract the type name from an alert method object as a string.

This function converts an alert method object (either an `AlertMethod` sum type
or one of its variant types) to its string representation and extracts just the
type name, removing any module qualifiers and parentheses.

# Arguments
- `alert_method`: An `AlertMethod` sum type or one of its variants
  (`DailyThreshold`, `MovingAverage`, `DailyThresholdMovingAverage`)

# Returns
- `String`: The unqualified type name (e.g., "DailyThreshold", "MovingAverage")

# Examples
```julia
getstring(DailyThreshold())  # Returns "DailyThreshold"
getstring(OutbreakDetectionCore.DailyThreshold())  # Returns "DailyThreshold"
```

# Implementation Notes
Uses a regex pattern `r"([^.(]+)\\("` to match the last component of a
potentially module-qualified type name before the opening parenthesis.

# See also
- [`AlertMethod`](@ref)
"""
alert_method_to_string(alert_method::AlertMethod) = match(
    r"([^.(]+)\(",
    string(
        LightSumTypes.variant(alert_method)
    )
).captures[1]

"""
    _construct_dirpath(alert_method, varargs...)

Internal dispatcher for directory path construction based on alert method type.

This function unwraps the `AlertMethod` sum type and dispatches to the
appropriate specialized path construction method.

# Arguments
- `alert_method`: An `AlertMethod` sum type
- `varargs...`: Variable arguments passed to the specialized methods
  (alertdirpath, testingdirpath, moving_average_lag)

# Returns
- `String`: Constructed directory path

# See also
- [`OutbreakDetectionSpecification`](@ref)
"""
function _construct_alert_method_dirpath(
        alert_method,
        varargs...
    )
    return _construct_alert_method_dirpath(
        LightSumTypes.variant(alert_method),
        varargs...
    )
end

"""
    _construct_dirpath(alert_method::DailyThreshold, alertdirpath,
                       testingdirpath, moving_average_lag)

Construct directory path for daily threshold alert method.

For the `DailyThreshold` method, the moving average lag is not used in the
path construction, so the path only includes alert and testing directories.

# Arguments
- `alert_method::DailyThreshold`: Alert method type indicator
- `alertdirpath`: Directory path component for alert parameters
- `testingdirpath`: Directory path component for testing parameters
- `moving_average_lag`: Moving average lag (unused for this method)

# Returns
- `String`: Joined path of alertdirpath and testingdirpath

# See also
- [`OutbreakDetectionSpecification`](@ref)
"""
function _construct_alert_method_dirpath(
        alert_method::DailyThreshold,
        alertdirpath,
        testingdirpath,
        moving_average_lag,
    )
    return joinpath(alertdirpath, testingdirpath)
end

"""
    _construct_dirpath(alert_method::Union{MovingAverage,
                       DailyThresholdMovingAverage}, alertdirpath,
                       testingdirpath, moving_average_lag)

Construct directory path for moving average-based alert methods.

For methods that use moving averages (`MovingAverage` and
`DailyThresholdMovingAverage`), the path includes a subdirectory for the
moving average lag parameter.

# Arguments
- `alert_method::Union{MovingAverage, DailyThresholdMovingAverage}`: Alert
  method type indicator
- `alertdirpath`: Directory path component for alert parameters
- `testingdirpath`: Directory path component for testing parameters
- `moving_average_lag`: Moving average lag value to include in path

# Returns
- `String`: Joined path including moving average lag subdirectory

# See also
- [`OutbreakDetectionSpecification`](@ref)
"""
function _construct_alert_method_dirpath(
        alert_method::Union{MovingAverage, DailyThresholdMovingAverage},
        alertdirpath,
        testingdirpath,
        moving_average_lag,
    )
    return joinpath(
        alertdirpath, "moveavglag_$(moving_average_lag)", testingdirpath
    )
end
