export calculate_positives, calculate_positives!,
    calculate_true_positives!, calculate_noise_positives!

"""
    calculate_positives(type_positive_function!, tested_vec, tlength, 
                       lag, testcharacteristic)

Calculate number of positive test results.
"""
function calculate_positives(
        type_positive_function!, tested_vec, tlength, lag, testcharacteristic
    )
    outvec = zeros(Int64, tlength)
    type_positive_function!(
        outvec, tested_vec, tlength, lag, testcharacteristic
    )
    return outvec
end

"""
    calculate_positives!(npos_vec, tested_vec, tlength, lag, tested_multiplier)

In-place calculation of positive test results with lag.
"""
function calculate_positives!(
        npos_vec, tested_vec, tlength, lag, tested_multiplier
    )
    @inbounds for day in eachindex(tested_vec)
        if day + lag <= tlength
            npos_vec[day + lag] = Int64(
                round(tested_vec[day] * tested_multiplier)
            )
        end
    end
    return nothing
end

"""
    calculate_noise_positives!(outvec, tested_vec, tlength, lag, spec)

Calculate false positive results from noise (1 - specificity).
"""
function calculate_noise_positives!(outvec, tested_vec, tlength, lag, spec)
    tested_multiplier = 1.0 - spec
    calculate_positives!(outvec, tested_vec, tlength, lag, tested_multiplier)
    return nothing
end

"""
    calculate_true_positives!(outvec, tested_vec, tlength, lag, sens)

Calculate true positive results (sensitivity).
"""
function calculate_true_positives!(outvec, tested_vec, tlength, lag, sens)
    calculate_positives!(outvec, tested_vec, tlength, lag, sens)
    return nothing
end
