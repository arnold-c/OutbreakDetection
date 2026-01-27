export match_outbreak_detection_bounds,
    get_matched_alert_indices,
    calculate_sensitivity,
    calculate_ppv

"""
    match_outbreak_detection_bounds(
        outbreakbounds::OutbreakThresholds,
        alertbounds::Thresholds,
        alert_filtering_strategy::AlertFilteringStrategy = AlertFilteringStrategy(AllAlerts())
    ) -> MatchedThresholds

Match outbreak periods to alert periods using memory-efficient sparse representation.

This function performs temporal matching between outbreak periods and alert periods,
returning a `MatchedThresholds` struct that uses sparse representation to minimize
memory usage.

# Matching Rule
Each alert is matched to **at most one outbreak** - specifically, the first outbreak
it overlaps with. This prevents double-counting alerts in PPV calculations and ensures
a clean one-to-one mapping.

# Arguments
- `outbreakbounds::OutbreakThresholds`: Outbreak periods with infection counts
- `alertbounds::Thresholds`: Alert periods from the detection system
- `alert_filtering_strategy::AlertFilteringStrategy`: Strategy for filtering alerts
  (default: `AllAlerts()`). Use `PostOutbreakStartAlerts()` to exclude alerts that
  start before the outbreak begins.

# Returns
- `MatchedThresholds`: Sparse representation containing:
  - Indices of outbreaks that have ≥1 alert
  - For each such outbreak, the indices of alerts that matched it
  - Total counts for metric calculation

# Matching Logic
An alert matches an outbreak if there is temporal overlap:
1. Alert starts within outbreak period (alert_lower >= outbreak_lower AND alert_lower <= outbreak_upper), OR
2. Alert starts before outbreak but extends into it (alert_lower <= outbreak_lower AND alert_upper > outbreak_lower)

If an alert overlaps multiple outbreaks, it is matched only to the first outbreak encountered.

When `PostOutbreakStartAlerts` filtering is used, only alerts where `alert_lower >= outbreak_lower`
are considered for matching (Case 1 only), ensuring all detection delays are non-negative.

# Example
```julia
outbreaks = OutbreakThresholds(
    lower_bounds = [10, 30, 50],
    upper_bounds = [20, 40, 60],
    duration = [11, 11, 11],
    num_infections_during_bounds = [100, 150, 200]
)

alerts = Thresholds(
    lower_bounds = [12, 35, 70],
    upper_bounds = [18, 38, 75],
    duration = [7, 4, 6]
)

matched = match_outbreak_detection_bounds(outbreaks, alerts)
# matched.outbreak_indices_with_alerts = [1, 2]  (outbreaks 1 and 2 detected, outbreak 3 missed)
# matched.alert_indices_per_outbreak = [[1], [2]]  (alert 1 → outbreak 1, alert 2 → outbreak 2)
# Alert 3 is a false alarm (no outbreak)

sensitivity = calculate_sensitivity(matched)  # 2/3 = 0.667
ppv = calculate_ppv(matched)                  # 2/3 = 0.667

# With PostOutbreakStartAlerts filtering (excludes pre-outbreak alerts)
matched_filtered = match_outbreak_detection_bounds(
    outbreaks, alerts, AlertFilteringStrategy(PostOutbreakStartAlerts())
)
```

# See Also
- [`MatchedThresholds`](@ref): Return type with sparse representation
- [`calculate_sensitivity`](@ref): Calculate proportion of outbreaks detected
- [`calculate_ppv`](@ref): Calculate proportion of alerts that are correct
- [`AlertFilteringStrategy`](@ref): Sum type for alert filtering strategies
"""
function match_outbreak_detection_bounds(
        outbreakbounds::OutbreakThresholds,
        alertbounds::Thresholds,
        alert_filtering_strategy::AlertFilteringStrategy = AlertFilteringStrategy(
            AllAlerts()
        ),
    )
    return _match_outbreak_detection_bounds(
        outbreakbounds,
        alertbounds,
        LightSumTypes.variant(alert_filtering_strategy),
    )
end

function _match_outbreak_detection_bounds(
        outbreakbounds::OutbreakThresholds,
        alertbounds::Thresholds,
        ::AllAlerts,
    )
    n_outbreaks = length(outbreakbounds.lower_bounds)
    n_alerts = length(alertbounds.lower_bounds)

    # Track which alerts have been matched (first outbreak only rule)
    matched_alerts = Set{Int64}()

    # For each outbreak, collect matching alert indices
    alert_indices_per_outbreak = [Int64[] for _ in eachindex(outbreakbounds.lower_bounds)]

    for outbreak_idx in eachindex(outbreakbounds.lower_bounds)
        outbreak_lower = outbreakbounds.lower_bounds[outbreak_idx]
        outbreak_upper = outbreakbounds.upper_bounds[outbreak_idx]

        for alert_idx in eachindex(alertbounds.lower_bounds)
            # Skip if alert already matched to a previous outbreak
            if alert_idx in matched_alerts
                continue
            end

            alert_lower = alertbounds.lower_bounds[alert_idx]
            alert_upper = alertbounds.upper_bounds[alert_idx]

            # Case 1: Alert starts within outbreak period
            if (alert_lower >= outbreak_lower && alert_lower <= outbreak_upper) ||
                    # Case 2: Alert starts before outbreak but extends into it
                    (alert_lower <= outbreak_lower && alert_upper > outbreak_lower)

                push!(alert_indices_per_outbreak[outbreak_idx], alert_idx)
                push!(matched_alerts, alert_idx)
            end
        end
    end

    # Filter to only outbreaks with at least one alert (sparse representation)
    outbreak_indices_with_alerts = findall(!isempty, alert_indices_per_outbreak)
    alert_indices_per_outbreak_filtered = alert_indices_per_outbreak[outbreak_indices_with_alerts]

    # Calculate matched counts for clarity (with n_outbreaks/alerts)
    # Also required for accuracy calculations so cache for performance
    n_matched_outbreaks = length(outbreak_indices_with_alerts)
    n_matched_alerts = length(matched_alerts)

    return MatchedThresholds(;
        outbreak_indices_with_alerts = outbreak_indices_with_alerts,
        alert_indices_per_outbreak = alert_indices_per_outbreak_filtered,
        n_matched_outbreaks = n_matched_outbreaks,
        n_matched_alerts = n_matched_alerts,
        n_outbreaks = n_outbreaks,
        n_alerts = n_alerts
    )
end

function _match_outbreak_detection_bounds(
        outbreakbounds::OutbreakThresholds,
        alertbounds::Thresholds,
        ::PostOutbreakStartAlerts,
    )
    n_outbreaks = length(outbreakbounds.lower_bounds)
    n_alerts = length(alertbounds.lower_bounds)

    # Track which alerts have been matched (first outbreak only rule)
    matched_alerts = Set{Int64}()

    # For each outbreak, collect matching alert indices
    alert_indices_per_outbreak = [Int64[] for _ in eachindex(outbreakbounds.lower_bounds)]

    for outbreak_idx in eachindex(outbreakbounds.lower_bounds)
        outbreak_lower = outbreakbounds.lower_bounds[outbreak_idx]
        outbreak_upper = outbreakbounds.upper_bounds[outbreak_idx]

        for alert_idx in eachindex(alertbounds.lower_bounds)
            # Skip if alert already matched to a previous outbreak
            if alert_idx in matched_alerts
                continue
            end

            alert_lower = alertbounds.lower_bounds[alert_idx]

            # PostOutbreakStartAlerts: Only match alerts that start on or after outbreak start
            # This ensures all detection delays are non-negative
            if alert_lower >= outbreak_lower && alert_lower <= outbreak_upper
                push!(alert_indices_per_outbreak[outbreak_idx], alert_idx)
                push!(matched_alerts, alert_idx)
            end
        end
    end

    # Filter to only outbreaks with at least one alert (sparse representation)
    outbreak_indices_with_alerts = findall(!isempty, alert_indices_per_outbreak)
    alert_indices_per_outbreak_filtered = alert_indices_per_outbreak[outbreak_indices_with_alerts]

    # Calculate matched counts for clarity (with n_outbreaks/alerts)
    # Also required for accuracy calculations so cache for performance
    n_matched_outbreaks = length(outbreak_indices_with_alerts)
    n_matched_alerts = length(matched_alerts)

    return MatchedThresholds(;
        outbreak_indices_with_alerts = outbreak_indices_with_alerts,
        alert_indices_per_outbreak = alert_indices_per_outbreak_filtered,
        n_matched_outbreaks = n_matched_outbreaks,
        n_matched_alerts = n_matched_alerts,
        n_outbreaks = n_outbreaks,
        n_alerts = n_alerts
    )
end

"""
    get_matched_alert_indices(matched::MatchedThresholds) -> Vector{Int64}

Get all alert indices that matched to at least one outbreak (flattened).

Returns a sorted vector of unique alert indices that were matched to outbreaks.
This is useful for calculating PPV or for further analysis of matched alerts.

# Arguments
- `matched::MatchedThresholds`: Matched thresholds structure

# Returns
- `Vector{Int64}`: Sorted vector of unique alert indices that matched outbreaks

# Example
```julia
matched = MatchedThresholds(
    outbreak_indices_with_alerts = [1, 3],
    alert_indices_per_outbreak = [[2, 5], [7, 9]],
    n_matched_outbreaks = 2,
    n_matched_alerts = 4,
    n_outbreaks = 5,
    n_alerts = 10
)

alert_indices = get_matched_alert_indices(matched)
# Returns: [2, 5, 7, 9]
```
"""
function get_matched_alert_indices(matched::MatchedThresholds)
    if isempty(matched.alert_indices_per_outbreak)
        return Int64[]
    end

    # Flatten and sort
    all_indices = reduce(vcat, matched.alert_indices_per_outbreak; init = Int64[])
    return sort(all_indices)
end
