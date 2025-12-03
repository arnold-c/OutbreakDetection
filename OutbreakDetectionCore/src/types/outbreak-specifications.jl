export OutbreakSpecification,
    OutbreakThresholdChars

"""
    OutbreakSpecification

Specification for outbreak detection thresholds.

# Fields
- `outbreak_threshold::Integer`: Threshold for outbreak detection
- `minimum_outbreak_duration::Integer`: Minimum duration to be considered outbreak
- `minimum_outbreak_size::Integer`: Minimum size to be considered outbreak
- `dirpath::AbstractString`: Directory path for output

# Constructor
    OutbreakSpecification(outbreak_threshold, minimum_outbreak_duration,
                         minimum_outbreak_size)

Creates an `OutbreakSpecification` with automatically generated directory path.
"""
struct OutbreakSpecification{T1 <: Integer, T2 <: AbstractString}
    outbreak_threshold::T1
    minimum_outbreak_duration::T1
    minimum_outbreak_size::T1
    dirpath::T2
end

function OutbreakSpecification(
        outbreak_threshold,
        minimum_outbreak_duration,
        minimum_outbreak_size
    )
    dirpath = joinpath(
        "min_outbreak_dur_$(minimum_outbreak_duration)",
        "min_outbreak_size_$(minimum_outbreak_size)",
        "outbreak_threshold_$(outbreak_threshold)",
    )

    return OutbreakSpecification(
        outbreak_threshold,
        minimum_outbreak_duration,
        minimum_outbreak_size,
        dirpath,
    )
end

"""
    OutbreakThresholdChars

Characteristics of outbreak detection performance for a given threshold.

Contains comprehensive metrics about detection accuracy, timing, and resource usage.
"""
struct OutbreakThresholdChars{
        T1 <: AbstractFloat,
        T2 <: Integer,
        T3 <: Vector{<:AbstractFloat},
        T4 <: Vector{<:Integer},
        T5 <: AbstractMatrix{<:Integer},
    }
    daily_sensitivity::T1
    daily_specificity::T1
    daily_ppv::T1
    daily_npv::T1
    accuracy::T1
    f1_score::T1
    matchedoutbreakbounds::T5
    noutbreaks::T2
    nalerts::T2
    outbreak_duration_vec::T4
    alert_duration_vec::T4
    detected_outbreak_size::T4
    missed_outbreak_size::T4
    n_true_outbreaks_detected::T2
    n_missed_outbreaks::T2
    n_correct_alerts::T2
    n_false_alerts::T2
    n_alerts_per_outbreak::T4
    period_sum_per_outbreak::T4
    perc_true_outbreaks_detected::T1
    perc_true_outbreaks_missed::T1
    falsealert_trueoutbreak_prop::T1
    correctalert_trueoutbreak_prop::T1
    trueoutbreak_alerts_prop::T1
    outbreaksmissed_alerts_prop::T1
    perc_alerts_false::T1
    perc_alerts_correct::T1
    detectiondelays::T4
    cases_before_alerts::T4
    cases_perc_before_alerts::T3
    cases_after_alerts::T4
    cases_perc_after_alerts::T3
    unavoidable_cases::T2
    avoidable_cases::T2
    n_outbreak_cases::T2
    n_tests::T2
    n_outbreak_tests::T2
    mean_noise_incidence_ratio::T1
    mean_poisson_noise::T1
    poisson_noise_prop::T1
    proportion_timeseries_in_outbreak::T1
    proportion_timeseries_in_alert::T1
    alert_outbreak_timeseries_prop_diff::T1
end
