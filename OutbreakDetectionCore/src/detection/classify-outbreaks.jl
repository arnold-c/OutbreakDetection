export classify_all_outbreaks,
    classify_outbreak

"""
    classify_all_outbreaks(
        incidence_vec,
        incidence_thresholds::Thresholds,
        minoutbreakdur,
        minoutbreaksize
    )

Classify threshold-crossing periods as outbreaks based on duration and size criteria.

This function takes periods where incidence exceeded a threshold and determines which
of these periods qualify as true outbreaks. It applies additional criteria for minimum
duration and total infection count to filter out brief spikes or small increases that
don't constitute meaningful outbreak events.

For each threshold-crossing period, the function:
1. Calculates the total number of infections during the period
2. Applies outbreak classification criteria using [`classify_outbreak`](@ref)
3. Retains only periods that meet both duration and size requirements
4. Returns outbreak-specific threshold information

# Arguments
- `incidence_vec`: Vector of daily incidence counts for a single simulation
- `incidence_thresholds::Thresholds`: Threshold periods with bounds and durations
- `minoutbreakdur`: Minimum consecutive days required for outbreak classification
- `minoutbreaksize`: Minimum total infections required for outbreak classification

# Returns
- `OutbreakThresholds`: Filtered threshold information containing only true outbreaks:
  - `lower_bounds`: Starting indices of outbreak periods
  - `upper_bounds`: Ending indices of outbreak periods
  - `duration`: Duration of each outbreak period
  - `num_infections_during_bounds`: Total infections during each outbreak period

# Examples
```julia
# Example threshold periods from incidence analysis
thresholds = Thresholds([10, 50], [15, 60], [6, 11])
incidence = rand(1:20, 100)

# Classify with outbreak criteria
outbreaks = classify_all_outbreaks(incidence, thresholds, 7, 50)

# Check which periods qualified as outbreaks
println("Outbreaks detected: ", length(outbreaks.lower_bounds))
```

# See Also
- [`classify_outbreak`](@ref): Function that applies classification criteria to individual periods
- [`Thresholds`](@ref): Input type containing threshold period information
- [`OutbreakThresholds`](@ref): Return type with outbreak-specific information
"""
function classify_all_outbreaks(
        incidence_vec,
        incidence_thresholds::Thresholds,
        minoutbreakdur,
        minoutbreaksize,
    )
    outbreak_status_vec = similar(incidence_thresholds.lower_bounds)
    ninf_during_bounds_vec = similar(incidence_thresholds.lower_bounds)

    for i in eachindex(incidence_thresholds.lower_bounds)
        @inbounds begin
            local lower_bound = incidence_thresholds.lower_bounds[i]
            local upper_bound = incidence_thresholds.upper_bounds[i]
            local duration = incidence_thresholds.duration[i]
            local ninf_during_bounds = sum(
                @view(incidence_vec[lower_bound:upper_bound])
            )
        end

        local outbreak_status = classify_outbreak(
            ninf_during_bounds,
            duration,
            minoutbreakdur,
            minoutbreaksize,
        )

        outbreak_status_vec[i] = outbreak_status
        ninf_during_bounds_vec[i] = ninf_during_bounds

    end

    outbreak_idx = findall(isequal(1), outbreak_status_vec)

    return OutbreakThresholds(
        Vector{Int64}(incidence_thresholds.lower_bounds[outbreak_idx]),
        Vector{Int64}(incidence_thresholds.upper_bounds[outbreak_idx]),
        Vector{Int64}(incidence_thresholds.duration[outbreak_idx]),
        Vector{Int64}(ninf_during_bounds_vec[outbreak_idx])
    )
end

"""
    classify_outbreak(
        ninf_during_bounds,
        duration,
        minoutbreakdur,
        minoutbreaksize
    )

Classify a single threshold period as an outbreak based on duration and size criteria.

This function applies the binary classification logic to determine whether a period
of above-threshold incidence qualifies as a true outbreak. A period is classified
as an outbreak if it meets both the minimum duration requirement (consecutive days)
and the minimum size requirement (total infections).

# Arguments
- `ninf_during_bounds`: Total number of infections during the threshold period
- `duration`: Number of consecutive days in the threshold period
- `minoutbreakdur`: Minimum consecutive days required for outbreak classification
- `minoutbreaksize`: Minimum total infections required for outbreak classification

# Returns
- `Int`: Binary classification result:
  - `1`: Period qualifies as an outbreak (meets both duration and size criteria)
  - `0`: Period does not qualify as an outbreak

# Examples
```julia
# Example: Period with 45 infections over 10 days
is_outbreak = classify_outbreak(45, 10, 7, 40)  # Returns 1 (outbreak)

# Example: Period with 30 infections over 5 days
is_outbreak = classify_outbreak(30, 5, 7, 40)   # Returns 0 (too short)

# Example: Period with 35 infections over 8 days
is_outbreak = classify_outbreak(35, 8, 7, 40)   # Returns 0 (too few infections)
```

# See Also
- [`classify_all_outbreaks`](@ref): Function that applies this to multiple periods
- [`OutbreakSpecification`](@ref): Type that defines the classification criteria
"""
function classify_outbreak(
        ninf_during_bounds,
        duration,
        minoutbreakdur,
        minoutbreaksize
    )
    if duration >= minoutbreakdur && ninf_during_bounds >= minoutbreaksize
        return 1
    end
    return 0
end
