export OutbreakDetectionSpecification, AlertMethod

# Define AlertMethod sum type variants
struct DailyThreshold end
struct MovingAverage end
struct DailyThresholdMovingAverage end

# Create the sum type
const AlertMethod = Union{DailyThreshold, MovingAverage, DailyThresholdMovingAverage}

"""
    AlertMethod

Method used for outbreak alert detection.

Available methods:
- `DailyThreshold`: Daily threshold only
- `MovingAverage`: Moving average only
- `DailyThresholdMovingAverage`: Both daily threshold and moving average
"""
AlertMethod

# String constructor for backwards compatibility
function alert_method_from_string(method_name::AbstractString)
    if method_name == "dailythreshold"
        return DailyThreshold()
    elseif method_name == "movingavg"
        return MovingAverage()
    elseif method_name == "dailythreshold_movingavg"
        return DailyThresholdMovingAverage()
    else
        available_methods = ["dailythreshold", "movingavg", "dailythreshold_movingavg"]
        error(
            "$(method_name) is not a valid alert method. It must be one of $(available_methods)"
        )
    end
end

# Display methods for printing
Base.show(io::IO, ::DailyThreshold) = print(io, "dailythreshold")
Base.show(io::IO, ::MovingAverage) = print(io, "movingavg")
Base.show(io::IO, ::DailyThresholdMovingAverage) =
    print(io, "dailythreshold_movingavg")

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
        TReal <: Real, T1 <: Integer, T2 <: AbstractFloat, T3 <: AlertMethod, T4 <: AbstractString,
    }
    alert_threshold::TReal
    moving_average_lag::T1
    percent_visit_clinic::T2
    percent_clinic_tested::T2
    percent_tested::T2
    alert_method::T3
    dirpath::T4
end

# Helper functions for directory path construction with dispatch
function _construct_dirpath(
        ::DailyThreshold,
        alertdirpath,
        testingdirpath,
        moving_average_lag,
    )
    return joinpath(alertdirpath, testingdirpath)
end

function _construct_dirpath(
        ::Union{MovingAverage, DailyThresholdMovingAverage},
        alertdirpath,
        testingdirpath,
        moving_average_lag,
    )
    return joinpath(
        alertdirpath, "moveavglag_$(moving_average_lag)", testingdirpath
    )
end

function OutbreakDetectionSpecification(
        alert_threshold,
        moving_average_lag,
        percent_visit_clinic,
        percent_clinic_tested,
        alert_method,
    )
    # Convert string to AlertMethod if needed
    alert_method_typed = alert_method isa AbstractString ? alert_method_from_string(alert_method) : alert_method

    alertdirpath = joinpath(
        "alertmethod_$(alert_method_typed)", "alertthreshold_$(alert_threshold)"
    )
    testingdirpath = joinpath(
        "perc_visit_clinic_$(percent_visit_clinic)",
        "perc_clinic_tested_$(percent_clinic_tested)",
    )

    dirpath = _construct_dirpath(
        alert_method_typed, alertdirpath, testingdirpath, moving_average_lag
    )

    return OutbreakDetectionSpecification(
        alert_threshold,
        moving_average_lag,
        percent_visit_clinic,
        percent_clinic_tested,
        percent_visit_clinic * percent_clinic_tested,
        alert_method_typed,
        dirpath,
    )
end
