export OutbreakDetectionSpecification, AlertMethod

"""
    AlertMethod

Method used for outbreak alert detection.

# Fields
- `method_name::AbstractString`: Name of the alert method

Available methods:
- "dailythreshold": Daily threshold only
- "movingavg": Moving average only
- "dailythreshold_movingavg": Both daily threshold and moving average
"""
struct AlertMethod{T1 <: AbstractString}
    method_name::T1
    function AlertMethod(method_name::T1) where {T1 <: AbstractString}
        available_test_methods = [
            "dailythreshold", "movingavg", "dailythreshold_movingavg",
        ]
        if !in(method_name, available_test_methods)
            error(
                "$(method_name) is not a valid test method. It must be one of $(available_test_methods)"
            )
        end
        return new{T1}(method_name)
    end
end

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

function OutbreakDetectionSpecification(
        alert_threshold,
        moving_average_lag,
        percent_visit_clinic,
        percent_clinic_tested,
        alert_method,
    )
    alertdirpath = joinpath(
        "alertmethod_$(alert_method)", "alertthreshold_$(alert_threshold)"
    )
    testingdirpath = joinpath(
        "perc_visit_clinic_$(percent_visit_clinic)",
        "perc_clinic_tested_$(percent_clinic_tested)",
    )

    dirpath = Match.@match alert_method begin
        "dailythreshold" => joinpath(
            alertdirpath,
            testingdirpath,
        )
        _ => joinpath(
            alertdirpath,
            "moveavglag_$(moving_average_lag)",
            testingdirpath,
        )
    end

    return OutbreakDetectionSpecification(
        alert_threshold,
        moving_average_lag,
        percent_visit_clinic,
        percent_clinic_tested,
        percent_visit_clinic * percent_clinic_tested,
        AlertMethod(alert_method),
        dirpath,
    )
end
