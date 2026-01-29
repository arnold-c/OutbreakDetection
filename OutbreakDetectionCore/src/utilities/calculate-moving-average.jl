export calculate_movingavg,
    calculate_movingavg!

function calculate_movingavg(
        vec_of_vecs::Vector{Vector{T}},
        avglag
    ) where {T <: Number}
    # Use Int64 with rounding to match old version behavior
    outvecs = Vector{Vector{Int64}}(undef, length(vec_of_vecs))
    for i in eachindex(vec_of_vecs)
        outvecs[i] = calculate_movingavg(
            vec_of_vecs[i],
            avglag
        )
    end
    return outvecs
end

"""
    calculate_movingavg(invec, avglag)

Calculate moving average of a vector.
Returns Int64 values (rounded) to match old version behavior.
"""
function calculate_movingavg(invec, avglag)
    # Use Int64 with rounding to match old version behavior
    outvec = zeros(Int64, size(invec, 1))

    calculate_movingavg!(outvec, invec, avglag)

    return outvec
end

"""
    calculate_movingavg!(outvec, invec, avglag)

In-place calculation of moving average.
"""
function calculate_movingavg!(outvec, invec, avglag)
    if avglag == 0
        outvec .= invec
        return nothing
    end

    @inbounds for day in eachindex(invec)
        outvec[day] = _calculate_int_daily_movingavg(invec, day, avglag)
    end
    return nothing
end

function _calculate_int_daily_movingavg(invec, day, avglag)
    @inline moveavg_daystart = _calculate_daily_movingavg_startday(day, avglag)
    return Int64(round(StatsBase.mean(@view(invec[moveavg_daystart:day]))))
end

function _calculate_daily_movingavg_startday(day, avglag)
    if day < avglag
        moveavg_daystart = 1
    else
        moveavg_daystart = day - avglag + 1
    end
    return moveavg_daystart
end
