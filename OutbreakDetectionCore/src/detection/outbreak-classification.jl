export classify_all_outbreaks!, filter_only_outbreaks,
    classify_outbreak, calculate_outbreak_duration,
    calculate_outbreak_size

"""
    classify_all_outbreaks!(outbreakstatus_vec, all_thresholds_arr, 
                            incidence_vec, minoutbreakdur, minoutbreaksize)

Classify all potential outbreaks based on duration and size criteria.
"""
function classify_all_outbreaks!(
        outbreakstatus_vec,
        all_thresholds_arr,
        incidence_vec,
        minoutbreakdur,
        minoutbreaksize,
    )
    for (row, (lower, upper)) in
        pairs(eachrow(@view(all_thresholds_arr[:, 1:2])))
        all_thresholds_arr[row, 3] = calculate_outbreak_duration(lower, upper)

        all_thresholds_arr[row, 4] = calculate_outbreak_size(
            incidence_vec, lower, upper
        )

        all_thresholds_arr[row, 5] = classify_outbreak(
            all_thresholds_arr[row, 3],
            minoutbreakdur,
            all_thresholds_arr[row, 4],
            minoutbreaksize,
        )

        @view(outbreakstatus_vec[lower:upper]) .= @view(
            all_thresholds_arr[row, 5]
        )
    end

    return nothing
end

"""
    classify_all_outbreaks!(outbreakstatus_vec, all_thresholds_arr, 
                            incidence_vec, outbreak_specification)

Classify all potential outbreaks using an OutbreakSpecification struct.
"""
function classify_all_outbreaks!(
        outbreakstatus_vec,
        all_thresholds_arr,
        incidence_vec,
        outbreak_specification::OutbreakSpecification,
    )
    return classify_all_outbreaks!(
        outbreakstatus_vec,
        all_thresholds_arr,
        incidence_vec,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )
end

"""
    calculate_outbreak_duration(lower_time, upper_time)

Calculate duration of an outbreak period.
"""
function calculate_outbreak_duration(outbreak_thresholds_row)
    return calculate_outbreak_duration(
        outbreak_thresholds_row[1], outbreak_thresholds_row[2]
    )
end

function calculate_outbreak_duration(lower_time, upper_time)
    return upper_time - lower_time + 1
end

function calculate_outbreak_duration!(outbreak_thresholds)
    outbreak_thresholds[:, 3] .=
        calculate_outbreak_duration.(
        outbreak_thresholds[:, 1], outbreak_thresholds[:, 2]
    )
    return nothing
end

"""
    calculate_outbreak_size(incidence_vec, lower_time, upper_time)

Calculate total size (cases) of an outbreak period.
"""
function calculate_outbreak_size(incidence_vec, lower_time, upper_time)
    return sum(@view(incidence_vec[lower_time:upper_time]))
end

"""
    filter_only_outbreaks(all_thresholds_arr)

Filter to keep only periods classified as true outbreaks.
"""
function filter_only_outbreaks(all_thresholds_arr)
    return @view(
        all_thresholds_arr[(all_thresholds_arr[:, 5] .== 1), 1:4]
    )
end

"""
    calculate_period_sum(incvec)

Calculate sum of incidence over a period.
"""
function calculate_period_sum(incvec)
    return sum(incvec)
end

"""
    classify_outbreak(outbreakdur, minoutbreakdur, periodsumvec, minoutbreaksize)

Classify whether a period qualifies as an outbreak (1) or not (0).
"""
function classify_outbreak(
        lower_time,
        upper_time,
        minoutbreakdur,
        periodsumvec,
        minoutbreaksize,
    )
    return classify_outbreak(
        calculate_outbreak_duration(lower_time, upper_time),
        minoutbreakdur,
        periodsumvec,
        minoutbreaksize,
    )
end

function classify_outbreak(
        outbreakdur, minoutbreakdur, periodsumvec, minoutbreaksize
    )
    if outbreakdur >= minoutbreakdur && periodsumvec >= minoutbreaksize
        return 1
    end
    return 0
end
