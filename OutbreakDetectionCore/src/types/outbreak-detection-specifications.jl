export OutbreakDetectionSpecification

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
    alert_method_string = alert_method_to_string(alert_method)
    alertdirpath = joinpath(
        "alertmethod_$(alert_method_string)", "alertthreshold_$(alert_threshold)"
    )
    testingdirpath = joinpath(
        "perc_visit_clinic_$(percent_visit_clinic)",
        "perc_clinic_tested_$(percent_clinic_tested)",
    )

    dirpath = _construct_alert_method_dirpath(
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
