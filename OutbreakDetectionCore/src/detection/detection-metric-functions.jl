export calculate_accuracy,
    calculate_sensitivity,
    calculate_ppv,
    arithmetic_mean,
    calculate_f_beta_score,
    calculate_detection_delay,
    calculate_unavoidable_cases,
    get_alert_duration,
    get_outbreak_duration,
    calculate_proportion_timeseries_in_alert,
    calculate_proportion_timeseries_in_outbreak

# Dispatch for accuracy metric
function calculate_accuracy(metric::AccuracyMetric, ppv, sensitivity)
    return _calculate_accuracy(LightSumTypes.variant(metric), ppv, sensitivity)
end

function _calculate_accuracy(::BalancedAccuracy, ppv, sensitivity)
    return arithmetic_mean(ppv, sensitivity)
end

function _calculate_accuracy(::F1, ppv, sensitivity)
    return calculate_f_beta_score(ppv, sensitivity)
end

"""
    arithmetic_mean(precision, recall)

Generic formula to calculate the arithmetic mean. Used for the `accuracy`
measure for outbreak detection.

  - Precision = PPV (% alerts that are correct)
  - Recall = sensitivity (% of outbreaks detected)

Implement as `arithmetic_mean(perc_alerts_correct, perc_true_outbreaks_detected)`
"""
function arithmetic_mean(precision, recall)
    return NaNMath.mean([precision, recall])
end

"""
    calculate_f_beta_score(precision, recall; beta = 1)

Generic formula to calculate the F-score. When beta=1, calculates the F1 score,
which weights precision and recall equally (the harmonic mean). beta=2 weights
precision more than recall, and beta=0.5 weights recall more.

  - Precision = PPV (% alerts that are correct)
  - Recall = sensitivity (% of outbreaks detected)

Implement as `calculate_f_beta_score(perc_alerts_correct, perc_true_outbreaks_detected; beta = 1)`
"""
function calculate_f_beta_score(
        precision, recall; beta = 1
    )
    f_beta =
        (1 + beta^2) * (precision * recall) / ((beta^2 * precision) + recall)

    if isnan(f_beta)
        f_beta = 0.0
    end

    return f_beta
end

"""
    calculate_sensitivity(matched::MatchedThresholds) -> Float64

Calculate sensitivity: proportion of outbreaks that have at least one matching alert.

Sensitivity measures the detection rate - what proportion of true outbreaks were
detected by the alert system. Also known as recall or true positive rate.

# Arguments
- `matched::MatchedThresholds`: Matched thresholds structure

# Returns
- `Float64`: Sensitivity in range [0, 1], or `NaN` if there are no outbreaks

# Formula
```
sensitivity = n_outbreaks_detected / n_total_outbreaks
```

# Example
```julia
matched = MatchedThresholds(
    outbreak_indices_with_alerts = [1, 3],  # 2 outbreaks detected
    alert_indices_per_outbreak = [[2], [7]],
    n_matched_outbreaks = 2,
    n_matched_alerts = 2,
    n_outbreaks = 5,  # out of 5 total
    n_alerts = 10
)

sensitivity = calculate_sensitivity(matched)
# Returns: 2/5 = 0.4
```

# See Also
- [`calculate_ppv`](@ref): Calculate positive predictive value
- [`MatchedThresholds`](@ref): Input structure
"""
function calculate_sensitivity(matched::MatchedThresholds)
    if matched.n_outbreaks == 0
        return NaN
    end

    return matched.n_matched_outbreaks / matched.n_outbreaks
end

"""
    calculate_ppv(matched::MatchedThresholds) -> Float64

Calculate PPV (Positive Predictive Value): proportion of alerts that matched an outbreak.

PPV measures the precision of the alert system - what proportion of alerts were
true positives (correctly identified an outbreak). Also known as precision.

# Arguments
- `matched::MatchedThresholds`: Matched thresholds structure

# Returns
- `Float64`: PPV in range [0, 1], or `NaN` if there are no alerts

# Formula
```
ppv = n_correct_alerts / n_total_alerts
```

# Example
```julia
matched = MatchedThresholds(
    outbreak_indices_with_alerts = [1, 3],
    alert_indices_per_outbreak = [[2, 5], [7]],  # 3 alerts matched
    n_matched_outbreaks = 2,
    n_matched_alerts = 3,
    n_outbreaks = 5,
    n_alerts = 10  # out of 10 total
)

ppv = calculate_ppv(matched)
# Returns: 3/10 = 0.3
```

# See Also
- [`calculate_sensitivity`](@ref): Calculate sensitivity (recall)
- [`MatchedThresholds`](@ref): Input structure
"""
function calculate_ppv(matched::MatchedThresholds)
    if matched.n_alerts == 0
        return NaN
    end

    return matched.n_matched_alerts / matched.n_alerts
end

"""
    calculate_detection_delay(matched::MatchedThresholds, outbreak_thresholds::OutbreakThresholds, alert_thresholds::Thresholds) -> Vector{Int64}

Calculate the time between the start of an outbreak and when the first matched alert occurs for each matched outbreak.

# Arguments
- `matched::MatchedThresholds`: Matched outbreak-alert pairs
- `outbreak_thresholds::OutbreakThresholds`: Outbreak bounds and metadata
- `alert_thresholds::Thresholds`: Alert bounds

# Returns
- `Vector{Int64}`: Detection delay for each matched outbreak in time units, or empty vector if no outbreaks were detected

# Formula
For each detected outbreak, the delay is:
```
delay = first_alert_start - outbreak_start
```
where `first_alert_start` is the lower bound of the first alert that matched the outbreak.

# Example
```julia
# Outbreak 1: days 10-20, first alert at day 15 → delay = 5
# Outbreak 3: days 30-40, first alert at day 32 → delay = 2
# Returns: [5, 2]
```
"""
function calculate_detection_delay(
        matched::MatchedThresholds,
        outbreak_thresholds::OutbreakThresholds,
        alert_thresholds::Thresholds,
    )
    if matched.n_matched_outbreaks == 0
        return Int64[]
    end

    delays = Vector{Int64}(undef, matched.n_matched_outbreaks)

    for (i, outbreak_idx) in enumerate(matched.outbreak_indices_with_alerts)
        # Get the start of the outbreak
        outbreak_start = outbreak_thresholds.lower_bounds[outbreak_idx]

        # Get the first alert that matched this outbreak
        first_alert_idx = first(matched.alert_indices_per_outbreak[i])
        first_alert_start = alert_thresholds.lower_bounds[first_alert_idx]

        # Calculate delay
        delays[i] = first_alert_start - outbreak_start
    end

    return delays
end

"""
    calculate_unavoidable_cases(matched::MatchedThresholds, outbreak_thresholds::OutbreakThresholds, alert_thresholds::Thresholds, incidence_vec::Vector{Int64}) -> Vector{Int64}

Calculate the number of cases that occur between the start of an outbreak and when the first matched alert occurs for each outbreak.

For outbreaks with matched alerts, unavoidable cases are those that occurred between the outbreak start
and the first alert. For outbreaks with NO matched alerts, ALL cases during the outbreak are unavoidable.

# Arguments
- `matched::MatchedThresholds`: Matched outbreak-alert pairs
- `outbreak_thresholds::OutbreakThresholds`: Outbreak bounds and metadata
- `alert_thresholds::Thresholds`: Alert bounds
- `incidence_vec::Vector{Int64}`: Time series of infection counts

# Returns
- `Vector{Int64}`: Number of unavoidable cases for each outbreak, or empty vector if no outbreaks exist

# Note
**MISSING CONTEXT**: This function requires `incidence_vec` which is not currently available in the calling context.
You will need to add this to the function signature in `calculate_optimal_results`.

# Formula
For matched outbreaks:
```
unavoidable_cases = sum(incidence_vec[outbreak_start:first_alert_start-1])
```

For unmatched outbreaks:
```
unavoidable_cases = sum(incidence_vec[outbreak_start:outbreak_end])
```

# Example
```julia
# Outbreak 1 (matched): days 10-20, first alert day 15
# incidence_vec[10:14] = [5, 8, 12, 15, 10] → unavoidable = 50

# Outbreak 2 (unmatched): days 30-35, no alerts
# incidence_vec[30:35] = [3, 5, 7, 6, 4, 2] → unavoidable = 27

# Returns: [50, 27]
```
"""
function calculate_unavoidable_cases(
        matched::MatchedThresholds,
        outbreak_thresholds::OutbreakThresholds,
        alert_thresholds::Thresholds,
        incidence_vec::Vector{Int64}
    )
    if matched.n_outbreaks == 0
        return Int64[]
    end

    unavoidable_cases = Vector{Int64}(undef, matched.n_outbreaks)

    # Create a set of matched outbreak indices for fast lookup
    matched_outbreak_set = Set(matched.outbreak_indices_with_alerts)

    for outbreak_idx in eachindex(unavoidable_cases)
        outbreak_start = outbreak_thresholds.lower_bounds[outbreak_idx]
        outbreak_end = outbreak_thresholds.upper_bounds[outbreak_idx]

        if outbreak_idx in matched_outbreak_set
            # This outbreak was detected - count cases before first alert
            # Find the position in the matched arrays
            matched_position = findfirst(
                isequal(outbreak_idx), matched.outbreak_indices_with_alerts
            )
            first_alert_idx = first(matched.alert_indices_per_outbreak[matched_position])
            first_alert_start = alert_thresholds.lower_bounds[first_alert_idx]

            # Sum infections from outbreak start to just before first alert
            if first_alert_start > outbreak_start
                unavoidable_cases[outbreak_idx] = sum(
                    @view(incidence_vec[outbreak_start:(first_alert_start - 1)])
                )
            else
                # Alert started at or before outbreak (no unavoidable cases)
                unavoidable_cases[outbreak_idx] = 0
            end
        else
            # This outbreak was never detected - ALL cases are unavoidable
            unavoidable_cases[outbreak_idx] = sum(
                @view(incidence_vec[outbreak_start:outbreak_end])
            )
        end
    end

    return unavoidable_cases
end

#     unavoidable_cases = Vector{Int64}(undef, matched.n_outbreaks)
#
#     # Create a set of matched outbreak indices for fast lookup
#     matched_outbreak_set = Set(matched.outbreak_indices_with_alerts)
#
#     for outbreak_idx in eachindex(unavoidable_cases)
#         outbreak_start = outbreak_thresholds.lower_bounds[outbreak_idx]
#         outbreak_end = outbreak_thresholds.upper_bounds[outbreak_idx]
#
#         if outbreak_idx in matched_outbreak_set
#             # This outbreak was detected - count cases before first alert
#             # Find the position in the matched arrays
#             matched_position = findfirst(
#                 isequal(outbreak_idx), matched.outbreak_indices_with_alerts
#             )
#             first_alert_idx = first(matched.alert_indices_per_outbreak[matched_position])
#             first_alert_start = alert_thresholds.lower_bounds[first_alert_idx]
#
#             # Sum infections from outbreak start to just before first alert
#             if first_alert_start > outbreak_start
#                 unavoidable_cases[outbreak_idx] = sum(
#                     @view(incidence_vec[outbreak_start:(first_alert_start - 1)])
#                 )
#             else
#                 # Alert started at or before outbreak (no unavoidable cases)
#                 unavoidable_cases[outbreak_idx] = 0
#             end
#         else
#             # This outbreak was never detected - ALL cases are unavoidable
#             unavoidable_cases[outbreak_idx] = sum(
#                 @view(incidence_vec[outbreak_start:outbreak_end])
#             )
#         end
#     end
#
#     return NaNMath.mean(unavoidable_cases)
# end

"""
    get_alert_duration(alert_thresholds::Thresholds) -> Vector{Int64}

get the duration of all alerts (not just matched alerts).

# Arguments
- `alert_thresholds::Thresholds`: Alert bounds containing duration information

# Returns
- `Vector{Int64}`: Duration of each alert in time units, or empty vector if there are no alerts

# Example
```julia
alert_thresholds = Thresholds(
    lower_bounds = [5, 15, 25],
    upper_bounds = [10, 20, 30],
    duration = [6, 6, 6]
)

durations = get_alert_duration(alert_thresholds)
# Returns: [6, 6, 6]
```
"""
function get_alert_duration(alert_thresholds::Thresholds)
    return alert_thresholds.duration
end

"""
    get_outbreak_duration(outbreak_thresholds::OutbreakThresholds) -> Vector{Int64}

Get the duration of all outbreaks (not just matched outbreaks).

# Arguments
- `outbreak_thresholds::OutbreakThresholds`: Outbreak bounds containing duration information

# Returns
- `Vector{Int64}`: Duration of each outbreak in time units, or empty vector if there are no outbreaks

# Example
```julia
outbreak_thresholds = OutbreakThresholds(
    lower_bounds = [10, 30, 50],
    upper_bounds = [20, 45, 65],
    duration = [11, 16, 16],
    num_infections_during_bounds = [100, 200, 150]
)

durations = get_outbreak_duration(outbreak_thresholds)
# Returns: [11, 16, 16]
```
"""
function get_outbreak_duration(outbreak_thresholds::OutbreakThresholds)
    return outbreak_thresholds.duration
end

"""
    calculate_proportion_timeseries_in_alert(alert_thresholds::Thresholds, timeseries_length::Int64) -> Float64

Calculate the proportion of the simulation timeseries that is in an alert status.

# Arguments
- `alert_thresholds::Thresholds`: Alert bounds containing duration information
- `timeseries_length::Int64`: Total length of the simulation timeseries

# Returns
- `Float64`: Proportion of timeseries in alert status in range [0, 1]

# Formula
```
proportion = total_alert_time / timeseries_length
```
where `total_alert_time = sum(alert_thresholds.duration)`

# Example
```julia
alert_thresholds = Thresholds(
    lower_bounds = [5, 15],
    upper_bounds = [10, 20],
    duration = [6, 6]
)
timeseries_length = 100

proportion = calculate_proportion_timeseries_in_alert(alert_thresholds, timeseries_length)
# Returns: (6 + 6) / 100 = 0.12
```
"""
function calculate_proportion_timeseries_in_alert(
        alert_thresholds::Thresholds,
        timeseries_length::Int64
    )
    if timeseries_length == 0
        return NaN
    end

    total_alert_time = sum(alert_thresholds.duration; init = 0)

    return total_alert_time / timeseries_length
end

"""
    calculate_proportion_timeseries_in_outbreak(outbreak_thresholds::OutbreakThresholds, timeseries_length::Int64) -> Float64

Calculate the proportion of the simulation timeseries that is in an outbreak status.

# Arguments
- `outbreak_thresholds::OutbreakThresholds`: Outbreak bounds containing duration information
- `timeseries_length::Int64`: Total length of the simulation timeseries

# Returns
- `Float64`: Proportion of timeseries in outbreak status in range [0, 1]

# Formula
```
proportion = total_outbreak_time / timeseries_length
```
where `total_outbreak_time = sum(outbreak_thresholds.duration)`

# Example
```julia
outbreak_thresholds = OutbreakThresholds(
    lower_bounds = [10, 30],
    upper_bounds = [20, 45],
    duration = [11, 16],
    num_infections_during_bounds = [100, 200]
)
timeseries_length = 100

proportion = calculate_proportion_timeseries_in_outbreak(outbreak_thresholds, timeseries_length)
# Returns: (11 + 16) / 100 = 0.27
```
"""
function calculate_proportion_timeseries_in_outbreak(
        outbreak_thresholds::OutbreakThresholds,
        timeseries_length::Int64
    )
    if timeseries_length == 0
        return NaN
    end

    total_outbreak_time = sum(outbreak_thresholds.duration; init = 0)

    return total_outbreak_time / timeseries_length
end
