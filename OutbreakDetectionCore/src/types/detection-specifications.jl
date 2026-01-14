export TestPositiveContainer,
    SingleAlertTestPositiveContainer,
    DualAlertTestPositiveContainer,
    OutbreakDetectionSpecification,
    AlertMethod,
    DailyThreshold,
    MovingAverage,
    DailyThresholdMovingAverage,
    getstring

abstract type AbstractTestPositiveContainer end

Base.@kwdef struct DualAlertTestPositiveContainer <: AbstractTestPositiveContainer
    test_positives::Vector{Vector{Int64}}
    movingavg_test_positives::Vector{Vector{Float64}}
end

Base.@kwdef struct SingleAlertTestPositiveContainer{T <: Union{Int64, Float64}}
    test_positive_results::Vector{Vector{T}}
end

LightSumTypes.@sumtype TestPositiveContainer(SingleAlertTestPositiveContainer, DualAlertTestPositiveContainer) <: AbstractTestPositiveContainer

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
getstring(alert_method::AlertMethod) = match(
    r"([^.(]+)\(",
    string(
        LightSumTypes.variant(alert_method)
    )
).captures[1]


# String constructor for backwards compatibility
"""
    OutbreakDetectionSpecification

Specification for outbreak detection parameters.

# Fields
- `alert_threshold::Real`: Threshold for triggering an alert
- `moving_average_lag::Integer`: Lag for moving average calculation
- `percent_visit_clinic::AbstractFloat`: Proportion visiting clinic
- `percent_clinic_tested::AbstractFloat`: Proportion of clinic visitors tested
- `percent_tested::AbstractFloat`: Overall proportion tested (product of above)
- `alert_method::AlertMethod`: Method for alert detection
- `dirpath::AbstractString`: Directory path for output

# Constructor
    OutbreakDetectionSpecification(alert_threshold, moving_average_lag,
                                   percent_visit_clinic, percent_clinic_tested,
                                   alert_method)

Creates an `OutbreakDetectionSpecification` with automatically generated directory path.
"""
struct OutbreakDetectionSpecification{
        TReal <: Real, T1 <: Integer, T2 <: AbstractFloat, T3 <: AbstractString,
    }
    alert_threshold::TReal
    moving_average_lag::T1
    percent_visit_clinic::T2
    percent_clinic_tested::T2
    percent_tested::T2
    alert_method::AlertMethod
    dirpath::T3
end

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
function _construct_dirpath(
        alert_method,
        varargs...
    )
    return _construct_dirpath(
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
function _construct_dirpath(
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
function _construct_dirpath(
        alert_method::Union{MovingAverage, DailyThresholdMovingAverage},
        alertdirpath,
        testingdirpath,
        moving_average_lag,
    )
    return joinpath(
        alertdirpath, "moveavglag_$(moving_average_lag)", testingdirpath
    )
end

"""
    OutbreakDetectionSpecification(alert_threshold, moving_average_lag,
                                   percent_visit_clinic, percent_clinic_tested,
                                   alert_method)

Construct an `OutbreakDetectionSpecification` with automatically generated
directory path.

This constructor creates a complete outbreak detection specification by
computing derived fields (overall testing percentage and directory path) from
the provided parameters.

# Arguments
- `alert_threshold`: Threshold value for triggering an outbreak alert
- `moving_average_lag`: Number of time steps for moving average calculation
- `percent_visit_clinic`: Proportion of infected individuals who visit a clinic
  (0.0 to 1.0)
- `percent_clinic_tested`: Proportion of clinic visitors who receive diagnostic
  testing (0.0 to 1.0)
- `alert_method`: Alert detection method (an `AlertMethod` sum type or variant)

# Returns
- `OutbreakDetectionSpecification`: Fully constructed specification with
  computed fields

# Computed Fields
- `percent_tested`: Automatically computed as `percent_visit_clinic *
  percent_clinic_tested`
- `dirpath`: Automatically generated hierarchical directory path based on all
  parameters

# Directory Path Structure
The generated `dirpath` follows this structure:
- For `DailyThreshold`:
  `"alertmethod_<method>/alertthreshold_<threshold>/perc_visit_clinic_<pvc>/perc_clinic_tested_<pct>"`
- For `MovingAverage` or `DailyThresholdMovingAverage`:
  `"alertmethod_<method>/alertthreshold_<threshold>/moveavglag_<lag>/perc_visit_clinic_<pvc>/perc_clinic_tested_<pct>"`

# Examples
```julia
spec = OutbreakDetectionSpecification(
    5.0,                              # alert_threshold
    7,                                # moving_average_lag
    0.8,                              # percent_visit_clinic
    0.9,                              # percent_clinic_tested
    AlertMethod(DailyThreshold())     # alert_method
)
# spec.percent_tested will be 0.72 (0.8 * 0.9)
```

# See also
- [`OutbreakDetectionSpecification`](@ref): The struct definition
- [`AlertMethod`](@ref)
- [`_construct_dirpath`](@ref)
"""
function OutbreakDetectionSpecification(
        alert_threshold,
        moving_average_lag,
        percent_visit_clinic,
        percent_clinic_tested,
        alert_method,
    )
    alert_method_string = getstring(alert_method)
    alertdirpath = joinpath(
        "alertmethod_$(alert_method_string)", "alertthreshold_$(alert_threshold)"
    )
    testingdirpath = joinpath(
        "perc_visit_clinic_$(percent_visit_clinic)",
        "perc_clinic_tested_$(percent_clinic_tested)",
    )

    dirpath = _construct_dirpath(
        alert_method,
        alertdirpath,
        testingdirpath,
        moving_average_lag
    )

    return OutbreakDetectionSpecification(
        alert_threshold,
        moving_average_lag,
        percent_visit_clinic,
        percent_clinic_tested,
        percent_visit_clinic * percent_clinic_tested,
        alert_method,
        dirpath,
    )
end
