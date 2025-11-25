export detectoutbreak, detectoutbreak!

"""
    _detect_outbreak!(
    	outbreakvec,
    	testarr_view,
    	test_movingavg_view,
    	alertthreshold,
    	alert_method::AlertMethod
    )

Internal dispatcher for outbreak detection based on alert method type.

This function unwraps the `AlertMethod` sum type and dispatches to the
appropriate specialized detection method.

# Arguments
- `outbreakvec`: Pre-allocated vector to store outbreak detection results (0/1)
- `testarr_view`: View of the daily test positivity or incidence array
- `test_movingavg_view`: View of the moving average array
- `alertthreshold`: Threshold value for triggering an outbreak alert
- `alert_method::AlertMethod`: Sum type specifying the detection method

# Returns
- `nothing` (modifies `outbreakvec` in-place)

# See also
- [`detectoutbreak!`](@ref)
- [`AlertMethod`](@ref)
"""
function _detect_outbreak!(
        outbreakvec,
        testarr_view,
        test_movingavg_view,
        alertthreshold,
        alert_method::AlertMethod,
    )
    return _detect_outbreak!(
        outbreakvec,
        testarr_view,
        test_movingavg_view,
        alertthreshold,
        LightSumTypes.variant(alert_method),
    )
end

"""
    _detect_outbreak!(
    	outbreakvec,
    	testarr_view,
    	test_movingavg_view,
    	alertthreshold,
    	alert_method::MovingAverage
    )

Detect outbreaks using only the moving average threshold method.

This dispatch variant is called when the alert method is `MovingAverage`,
which compares only the moving average values against the threshold.

# Arguments
- `outbreakvec`: Pre-allocated vector to store outbreak detection results (0/1)
- `testarr_view`: View of the daily test positivity or incidence array (unused)
- `test_movingavg_view`: View of the moving average array
- `alertthreshold`: Threshold value for triggering an outbreak alert
- `alert_method::MovingAverage`: Alert method type indicator

# Returns
- `nothing` (modifies `outbreakvec` in-place)
"""
function _detect_outbreak!(
        outbreakvec,
        testarr_view,
        test_movingavg_view,
        alertthreshold,
        alert_method::MovingAverage,
    )
    return detectoutbreak!(outbreakvec, test_movingavg_view, alertthreshold)
end

"""
    _detect_outbreak!(
    	outbreakvec,
    	testarr_view,
    	test_movingavg_view,
    	alertthreshold,
    	alert_method::DailyThreshold
    )

Detect outbreaks using only the daily threshold method.

This dispatch variant is called when the alert method is `DailyThreshold`,
which compares only the daily values against the threshold.

# Arguments
- `outbreakvec`: Pre-allocated vector to store outbreak detection results (0/1)
- `testarr_view`: View of the daily test positivity or incidence array
- `test_movingavg_view`: View of the moving average array (unused)
- `alertthreshold`: Threshold value for triggering an outbreak alert
- `alert_method::DailyThreshold`: Alert method type indicator

# Returns
- `nothing` (modifies `outbreakvec` in-place)
"""
function _detect_outbreak!(
        outbreakvec,
        testarr_view,
        test_movingavg_view,
        alertthreshold,
        alert_method::DailyThreshold,
    )
    return detectoutbreak!(outbreakvec, testarr_view, alertthreshold)
end


"""
    _detect_outbreak!(
    	outbreakvec,
    	testarr_view,
    	test_movingavg_view,
    	alertthreshold,
    	alert_method::DailyThresholdMovingAverage
    )

Detect outbreaks using both daily and moving average threshold methods.

This dispatch variant is called when the alert method is
`DailyThresholdMovingAverage`, which triggers an alert if either the daily
value or the moving average exceeds the threshold.

# Arguments
- `outbreakvec`: Pre-allocated vector to store outbreak detection results (0/1)
- `testarr_view`: View of the daily test positivity or incidence array
- `test_movingavg_view`: View of the moving average array
- `alertthreshold`: Threshold value for triggering an outbreak alert
- `alert_method::DailyThresholdMovingAverage`: Alert method type indicator

# Returns
- `nothing` (modifies `outbreakvec` in-place)
"""
function _detect_outbreak!(
        outbreakvec,
        testarr_view,
        test_movingavg_view,
        alertthreshold,
        alert_method::DailyThresholdMovingAverage,
    )
    return detectoutbreak!(
        outbreakvec, testarr_view, test_movingavg_view, alertthreshold
    )
end

"""
    detectoutbreak(incvec, avgvec, threshold)

Detect outbreak based on daily and moving average thresholds.

Creates a new outbreak detection vector and populates it by comparing both
the daily incidence and moving average values against the threshold. An
outbreak is detected (value = 1) if either the daily value or moving average
exceeds the threshold.

# Arguments
- `incvec`: Vector of daily incidence values
- `avgvec`: Vector of moving average values
- `threshold`: Threshold value for triggering an outbreak alert

# Returns
- `outbreak`: Vector of Int64 values (0 or 1) indicating outbreak status at
  each time point

# See also
- [`detectoutbreak!`](@ref): In-place version of this function
"""
function detectoutbreak(incvec, avgvec, threshold)
    outbreak = zeros(Int64, length(incvec))

    detectoutbreak!(outbreak, incvec, avgvec, threshold)

    return outbreak
end

"""
    detectoutbreak!(outbreakvec, incvec, avgvec, threshold)

In-place outbreak detection using both daily and moving average thresholds.

Determines whether an outbreak has been detected at each time point by
comparing both the daily incidence and moving average values against the
threshold. An outbreak is detected (value = 1) if either condition is met.

# Arguments
- `outbreakvec`: Pre-allocated vector to store outbreak detection results (0/1)
- `incvec`: Vector of daily incidence values
- `avgvec`: Vector of moving average values of daily incidence
- `threshold`: Threshold value for triggering an outbreak alert

# Returns
- `nothing` (modifies `outbreakvec` in-place)

# Notes
The function uses element-wise operations for efficiency. At each time point,
`outbreakvec[i] = 1` if `incvec[i] >= threshold` OR `avgvec[i] >= threshold`,
otherwise `outbreakvec[i] = 0`.

# See also
- [`detectoutbreak`](@ref): Allocating version of this function
"""
function detectoutbreak!(outbreakvec, incvec, avgvec, threshold)
    @. outbreakvec = ifelse(incvec >= threshold || avgvec >= threshold, 1, 0)

    return nothing
end

"""
    detectoutbreak(infectionsvec, threshold)

Detect outbreak based on a single threshold.

Creates a new outbreak detection vector and populates it by comparing the
infection values against the threshold. An outbreak is detected (value = 1)
when the infection value meets or exceeds the threshold.

# Arguments
- `infectionsvec`: Vector of infection values (either daily incidence or
  moving average)
- `threshold`: Threshold value for triggering an outbreak alert

# Returns
- `outbreak`: Vector of Int64 values (0 or 1) indicating outbreak status at
  each time point

# Notes
This function can be used with either daily incidence values or moving
average values, depending on the desired detection method.

# See also
- [`detectoutbreak!`](@ref): In-place version of this function
"""
function detectoutbreak(infectionsvec, threshold)
    outbreak = zeros(Int64, length(infectionsvec))

    detectoutbreak!(outbreak, infectionsvec, threshold)

    return outbreak
end

"""
    detectoutbreak!(outbreakvec, infectionsvec, threshold)

In-place outbreak detection using a single threshold.

Determines whether an outbreak has been detected at each time point by
comparing the infection values against the threshold. An outbreak is
detected (value = 1) when the infection value meets or exceeds the threshold.

# Arguments
- `outbreakvec`: Pre-allocated vector to store outbreak detection results (0/1)
- `infectionsvec`: Vector of infection values (either daily incidence or
  moving average)
- `threshold`: Threshold value for triggering an outbreak alert

# Returns
- `nothing` (modifies `outbreakvec` in-place)

# Notes
The `infectionsvec` can be either:
- A vector of daily infections (for daily threshold method)
- A vector of moving average values (for moving average method)

The function uses element-wise operations for efficiency. At each time point,
`outbreakvec[i] = 1` if `infectionsvec[i] >= threshold`, otherwise
`outbreakvec[i] = 0`.

# See also
- [`detectoutbreak`](@ref): Allocating version of this function
"""
function detectoutbreak!(outbreakvec, infectionsvec, threshold)
    @. outbreakvec = ifelse(infectionsvec >= threshold, 1, 0)

    return nothing
end
