export Thresholds,
    OutbreakThresholds,
    MatchedThresholds

abstract type AbstractThresholds end

"""
    Thresholds <: AbstractThresholds

Basic threshold specification for epidemiological surveillance analysis.

This struct defines threshold parameters used in outbreak detection and effective reproduction
number (Reff) analysis, specifying the bounds and duration criteria for identifying significant
epidemiological events. It provides the fundamental threshold configuration for time-series
analysis of disease surveillance data.

The thresholds can represent:
- Outbreak periods: Time intervals when incidence exceeds a specified threshold
- Reff periods: Time intervals when the effective reproduction number (Reff) is >= 1

# Fields
- `lower_bounds::Vector{Int64}`: Starting time indices of threshold exceedance periods
- `upper_bounds::Vector{Int64}`: Ending time indices of threshold exceedance periods
- `duration::Vector{Int64}`: Duration (in time periods) of each threshold exceedance

# Constructor
    Thresholds(; lower_bounds, upper_bounds, duration)

# Example
```julia
# Create threshold specification for outbreak detection
outbreak_thresholds = Thresholds(
    lower_bounds = [5, 10, 15],
    upper_bounds = [20, 30, 40],
    duration = [16, 21, 26]  # upper - lower + 1
)

# Create threshold specification for Reff >= 1 periods
reff_thresholds = Reff_ge_than_one(Reff_vec)

# Access threshold parameters
println("Lower bounds: \$(outbreak_thresholds.lower_bounds)")
println("Duration requirements: \$(outbreak_thresholds.duration)")
```

# See Also
- [`AbstractThresholds`](@ref): Parent abstract type
- [`OutbreakThresholds`](@ref): Extended threshold specification with infection counts
- [`Reff_ge_than_one`](@ref): Function for calculating Reff >= 1 thresholds
"""
Base.@kwdef struct Thresholds <: AbstractThresholds
    lower_bounds::Vector{Int64}
    upper_bounds::Vector{Int64}
    duration::Vector{Int64}
end


"""
    OutbreakThresholds <: AbstractThresholds

Extended threshold specification for outbreak detection with infection count tracking.

This struct extends the basic threshold specification to include tracking of infection
counts during threshold exceedance periods. It provides comprehensive threshold
configuration for detailed outbreak analysis that requires monitoring both detection
criteria and the epidemiological burden during outbreak periods.

Tracking the infection counts allow for the classification of
whether a threshold exceedance qualifies as a true outbreak based on minimum
duration and size criteria.

# Fields
- `lower_bounds::Vector{Int64}`: Starting time indices of outbreak periods
- `upper_bounds::Vector{Int64}`: Ending time indices of outbreak periods
- `duration::Vector{Int64}`: Duration (in time periods) of each outbreak
- `num_infections_during_bounds::Vector{Int64}`: Total number of infections during each outbreak period

# Constructor
    OutbreakThresholds(; lower_bounds, upper_bounds, duration, num_infections_during_bounds)

# Example
```julia
# Create extended threshold specification with infection tracking
outbreak_thresholds = OutbreakThresholds(
    lower_bounds = [5, 10, 15],
    upper_bounds = [20, 30, 40],
    duration = [16, 21, 26],
    num_infections_during_bounds = [50, 150, 300]
)

# Access infection count data
println("Infections during bounds: \$(outbreak_thresholds.num_infections_during_bounds)")
println("Total outbreaks detected: \$(length(outbreak_thresholds.lower_bounds))")
```

# See Also
- [`AbstractThresholds`](@ref): Parent abstract type
- [`Thresholds`](@ref): Basic threshold specification without infection tracking
- [`OutbreakSpecification`](@ref): Outbreak definition parameters
- [`calculate_outbreak_thresholds`](@ref): Function for calculating outbreak thresholds
"""
Base.@kwdef struct OutbreakThresholds <: AbstractThresholds
    lower_bounds::Vector{Int64}
    upper_bounds::Vector{Int64}
    duration::Vector{Int64}
    num_infections_during_bounds::Vector{Int64}
end


"""
    MatchedThresholds <: AbstractThresholds

Memory-efficient sparse representation of matched alert-outbreak pairs.

This struct stores the results of matching alerts to outbreaks, using a sparse
representation that only stores outbreaks that have at least one matching alert.

# Matching Rule
An alert can match multiple outbreaks if it spans across them. Specifically, an alert
that was matched to outbreak N via Case 1 (alert starts within outbreak) can also match
outbreak N+1 via Case 2 (alert starts before outbreak but extends into it) if the alert
extends past the end of outbreak N into outbreak N+1.

# Fields
- `outbreak_indices_with_alerts::Vector{Int64}`: Indices of outbreaks that have ≥1 matching alert
- `alert_indices_per_outbreak::Vector{Vector{Int64}}`: For each outbreak in `outbreak_indices_with_alerts`,
  the indices of alerts that matched it (into the original alert bounds)
- `n_matched_outbreaks::Int64`: Number of outbreaks that have ≥1 matching alert (cached for clarity)
- `n_matched_alerts::Int64`: Total number of alerts that matched to an outbreak (cached for clarity)
- `n_outbreaks::Int64`: Total number of outbreaks (including those with no alerts)
- `n_alerts::Int64`: Total number of alerts (including those that didn't match any outbreak)

# Invariants
The struct validates that:
- `length(outbreak_indices_with_alerts) == length(alert_indices_per_outbreak)`
- `n_matched_outbreaks == length(outbreak_indices_with_alerts)`
- `n_matched_alerts == sum(length, alert_indices_per_outbreak)`

# Constructor
    MatchedThresholds(; outbreak_indices_with_alerts, alert_indices_per_outbreak, n_matched_outbreaks, n_matched_alerts, n_outbreaks, n_alerts)

# Example
```julia
# Scenario: 5 outbreaks, 10 alerts, only outbreaks 1 and 3 detected
matched = MatchedThresholds(
    outbreak_indices_with_alerts = [1, 3],
    alert_indices_per_outbreak = [[2, 5], [7]],  # outbreak 1 matched alerts 2,5; outbreak 3 matched alert 7
    n_matched_outbreaks = 2,  # 2 outbreaks detected
    n_matched_alerts = 3,     # 3 alerts matched (2 + 1)
    n_outbreaks = 5,
    n_alerts = 10
)

# Calculate metrics (now O(1) operations using cached counts)
sensitivity = calculate_sensitivity(matched)  # 2/5 = 0.4 (uses n_matched_outbreaks)
ppv = calculate_ppv(matched)                  # 3/10 = 0.3 (uses n_matched_alerts)
```

# Memory Efficiency
In low-sensitivity scenarios, this representation is significantly more memory-efficient
than storing full parallel vectors:
- 5 detected outbreaks out of 100: ~75% memory savings vs parallel vectors
- 50 detected outbreaks out of 100: ~20% memory savings vs parallel vectors

# See Also
- [`AbstractThresholds`](@ref): Parent abstract type
- [`Thresholds`](@ref): Basic threshold specification
- [`OutbreakThresholds`](@ref): Outbreak threshold specification with infection counts
- [`match_outbreak_detection_bounds`](@ref): Function that creates MatchedThresholds
- [`calculate_sensitivity`](@ref): Calculate proportion of outbreaks detected
- [`calculate_ppv`](@ref): Calculate proportion of alerts that are correct
"""
Base.@kwdef struct MatchedThresholds <: AbstractThresholds
    outbreak_indices_with_alerts::Vector{Int64}
    alert_indices_per_outbreak::Vector{Vector{Int64}}
    n_matched_outbreaks::Int64
    n_matched_alerts::Int64
    n_outbreaks::Int64
    n_alerts::Int64

    function MatchedThresholds(
            outbreak_indices_with_alerts::Vector{Int64},
            alert_indices_per_outbreak::Vector{Vector{Int64}},
            n_matched_outbreaks::Int64,
            n_matched_alerts::Int64,
            n_outbreaks::Int64,
            n_alerts::Int64
        )
        @assert length(outbreak_indices_with_alerts) == length(alert_indices_per_outbreak) "Length mismatch: outbreak_indices_with_alerts has $(length(outbreak_indices_with_alerts)) elements but alert_indices_per_outbreak has $(length(alert_indices_per_outbreak)) elements"
        @assert n_outbreaks >= 0 "n_outbreaks must be non-negative, got $n_outbreaks"
        @assert n_alerts >= 0 "n_alerts must be non-negative, got $n_alerts"
        @assert n_matched_outbreaks >= 0 "n_matched_outbreaks must be non-negative, got $n_matched_outbreaks"
        @assert n_matched_alerts >= 0 "n_matched_alerts must be non-negative, got $n_matched_alerts"
        @assert n_matched_outbreaks == length(outbreak_indices_with_alerts) "n_matched_outbreaks ($n_matched_outbreaks) must equal length(outbreak_indices_with_alerts) ($(length(outbreak_indices_with_alerts)))"
        @assert n_matched_alerts == sum(length, alert_indices_per_outbreak; init = 0) "n_matched_alerts ($n_matched_alerts) must equal sum of alert counts ($(sum(length, alert_indices_per_outbreak; init = 0)))"
        @assert all(
            idx -> 1 <= idx <= n_outbreaks, outbreak_indices_with_alerts
        ) "All outbreak indices must be in range [1, $n_outbreaks]"
        @assert all(
            alerts -> all(idx -> 1 <= idx <= n_alerts, alerts),
            alert_indices_per_outbreak
        ) "All alert indices must be in range [1, $n_alerts]"
        @assert allunique(outbreak_indices_with_alerts) "There are non-unique outbreak indices"
        # NOTE: Alert indices can be non-unique because an alert can match multiple outbreaks
        # if it spans across them (via Case 2: alert starts before outbreak but extends into it)

        return new(
            outbreak_indices_with_alerts,
            alert_indices_per_outbreak,
            n_matched_outbreaks,
            n_matched_alerts,
            n_outbreaks,
            n_alerts
        )
    end
end
