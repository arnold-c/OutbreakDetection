export calculate_test_positivity

"""
    calculate_test_positivity(true_positive_vec, total_test_vec, 
                              alert_vec, agg_days)

Calculate test positivity rate over aggregated time periods.

Returns a matrix with test positivity rates and outbreak status.
"""
function calculate_test_positivity(
        true_positive_vec, total_test_vec, alert_vec, agg_days
    )
    @views outvec = zeros(Float64, length(true_positive_vec) ÷ agg_days, 2)
    @inbounds for i in axes(outvec, 1)
        start_ind = 1 + (i - 1) * agg_days
        end_ind = start_ind + (agg_days - 1)

        @views total_test_sum = sum(total_test_vec[start_ind:end_ind])
        @views true_positive_sum = sum(true_positive_vec[start_ind:end_ind])
        @views num_outbreak_days = sum(alert_vec[start_ind:end_ind])
        agg_outbreak_status = num_outbreak_days >= agg_days / 2 ? 1 : 0

        outvec[i, 1] = true_positive_sum / total_test_sum
        outvec[i, 2] = agg_outbreak_status
    end
    return outvec
end
