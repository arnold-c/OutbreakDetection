export detectoutbreak, detectoutbreak!

"""
    detectoutbreak(incvec, avgvec, threshold)

Detect outbreak based on daily and moving average thresholds.
"""
function detectoutbreak(incvec, avgvec, threshold)
    outbreak = zeros(Int64, length(incvec))

    detectoutbreak!(outbreak, incvec, avgvec, threshold)

    return outbreak
end

"""
    detectoutbreak!(outbreakvec, incvec, avgvec, threshold)

Determines whether an outbreak has been detected based on an infection daily 
and moving average timeseries and a threshold.

`incvec` should be the daily incidence timeseries, and `avgvec` should be 
the moving average of the daily incidence timeseries.
"""
function detectoutbreak!(outbreakvec, incvec, avgvec, threshold)
    @. outbreakvec = ifelse(incvec >= threshold || avgvec >= threshold, 1, 0)

    return nothing
end

"""
    detectoutbreak(infectionsvec, threshold)

Detect outbreak based on a single threshold.
"""
function detectoutbreak(infectionsvec, threshold)
    outbreak = zeros(Int64, length(infectionsvec))

    detectoutbreak!(outbreak, infectionsvec, threshold)

    return outbreak
end

"""
    detectoutbreak!(outbreakvec, infectionsvec, threshold)

Determines whether an outbreak has been detected based on an infection 
timeseries and a threshold.

The `infectionsvec` can either be a vector of daily infections or a vector 
of the moving average of daily infections.
"""
function detectoutbreak!(outbreakvec, infectionsvec, threshold)
    @. outbreakvec = ifelse(infectionsvec >= threshold, 1, 0)

    return nothing
end
