# module DiagTestingFunctions
#
# export create_testing_arr, create_testing_arr!, calculate_tested!,
#     calculate_pos, calculate_pos!, calculate_movingavg, calculate_movingavg!,
#     detectoutbreak, detectoutbreak!, calculate_ot_characterstics,
#     calculate_noutbreaks, calculate_OutbreakThresholdChars,
#     run_OutbreakThresholdChars_creation, OutbreakThresholdChars_creation

using DrWatson
using StatsBase
using FreqTables
using ThreadsX
using FLoops
using StructArrays
using NaNMath: NaNMath

# include("detection-thresholds.jl")
# # using .DetectionThresholds
#
# include("structs.jl")
# using .ODStructs

function create_testing_arrs(
    incarr,
    noisearr,
    outbreak_detect_spec::OutbreakDetectionSpecification,
    individual_test_spec::IndividualTestSpecification,
)
    testarr = zeros(Int64, size(incarr, 1), 8, size(incarr, 3))
    testpos_vec = Vector{TestPositivity}(undef, size(incarr, 3))
    ntested_worker_vec = Vector{Int64}(undef, size(incarr, 1))

    create_testing_arrs!(
        testarr,
        testpos_vec,
        ntested_worker_vec,
        incarr,
        noisearr,
        outbreak_detect_spec.alert_threshold,
        outbreak_detect_spec.moving_average_lag,
        outbreak_detect_spec.percent_tested,
        individual_test_spec.test_result_lag,
        individual_test_spec.sensitivity,
        individual_test_spec.specificity,
    )

    return testarr, StructArray(testpos_vec)
end

function create_testing_arrs!(
    testarr,
    testpos_vec,
    ntested_worker_vec,
    incarr,
    noisearr,
    alertthreshold,
    moveavglag,
    perc_tested,
    testlag,
    testsens,
    testspec,
)
    tlength = size(testarr, 1)

    for sim in axes(incarr, 3)
        # Number of infectious individuals tested
        calculate_tested!(
            @view(testarr[:, 1, sim]), @view(incarr[:, 1, sim]), perc_tested
        )

        # Number of noise individuals tested
        calculate_tested!(
            @view(testarr[:, 2, sim]), @view(noisearr[:, 1, sim]), perc_tested
        )

        # Number of TOTAL individuals tested
        @turbo @. @views ntested_worker_vec .=
            testarr[:, 1, sim] + testarr[:, 2, sim]

        # Number of test positive INFECTED individuals
        calculate_true_positives!(
            @view(testarr[:, 3, sim]),
            @view(testarr[:, 1, sim]),
            tlength,
            testlag,
            testsens,
        )

        # Number of test positive NOISE individuals
        calculate_noise_positives!(
            @view(testarr[:, 4, sim]),
            @view(testarr[:, 2, sim]),
            tlength,
            testlag,
            testspec,
        )

        # Number of test positive TOTAL individuals
        @. testarr[:, 5, sim] =
            @view(testarr[:, 3, sim]) + @view(testarr[:, 4, sim])

        # Calculate moving average of TOTAL test positives
        calculate_movingavg!(
            @view(testarr[:, 6, sim]),
            @view(testarr[:, 5, sim]),
            moveavglag
        )

        # TOTAL Test positive individuals trigger outbreak response
        detectoutbreak!(
            @view(testarr[:, 7, sim]),
            @view(testarr[:, 5, sim]),
            @view(testarr[:, 6, sim]),
            alertthreshold,
        )

        # Triggered outbreak equal to actual outbreak status
        @. testarr[:, 8, sim] =
            @view(testarr[:, 7, sim]) == @view(incarr[:, 3, sim])

        # Posterior prob of infectious / total test tests performed
        testpos_vec[sim] = TestPositivity(
            @view(testarr[:, 5, sim]),
            ntested_worker_vec,
            @view(testarr[:, 7, sim])
        )
    end

    return nothing
end

function calculate_tested!(outvec, invec, perc_tested)
    @. outvec = round(invec * perc_tested)
end

function calculate_noise_positives!(outvec, tested_vec, tlength, lag, spec)
    tested_multiplier = 1.0 - spec
    calculate_positives!(outvec, tested_vec, tlength, lag, tested_multiplier)
    return nothing
end

function calculate_true_positives!(outvec, tested_vec, tlength, lag, sens)
    calculate_positives!(outvec, tested_vec, tlength, lag, sens)
    return nothing
end

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

function calculate_movingavg(invec, avglag)
    outvec = zeros(Float64, size(invec, 1))

    calculate_movingavg!(outvec, invec, avglag)

    return outvec
end

function calculate_movingavg(
    invec::T1, avglag
) where {T1<:AbstractArray{Integer}}
    outvec = zeros(eltype(invec), size(invec, 1))

    calculate_movingavg!(outvec, invec, avglag)

    return outvec
end

function calculate_movingavg!(outvec, invec, avglag)
    if avglag == 0
        outvec .= invec
        return nothing
    end

    @inbounds for day in eachindex(invec)
        outvec[day] = calculate_float_daily_movingavg(invec, day, avglag)
    end
    return nothing
end

function calculate_movingavg!(
    outvec::T1, invec, avglag
) where {T1<:AbstractArray{<:Integer}}
    if avglag == 0
        outvec .= invec
        return nothing
    end

    @inbounds for day in eachindex(invec)
        outvec[day] = calculate_int_daily_movingavg(invec, day, avglag)
    end
    return nothing
end

function calculate_float_daily_movingavg(invec, day, avglag)
    @inline moveavg_daystart = calculate_daily_movingavg_startday(day, avglag)
    return mean(@view(invec[moveavg_daystart:day]))
end

function calculate_int_daily_movingavg(invec, day, avglag)
    @inline moveavg_daystart = calculate_daily_movingavg_startday(day, avglag)
    return Int64(round(mean(@view(invec[moveavg_daystart:day]))))
end

function calculate_daily_movingavg_startday(day, avglag)
    if day < avglag
        moveavg_daystart = 1
    else
        moveavg_daystart = day - avglag + 1
    end
    return moveavg_daystart
end

function detectoutbreak(incvec, avgvec, threshold)
    outbreak = zeros(Int64, length(incvec))

    detectoutbreak!(outbreak, incvec, avgvec, threshold)

    return outbreak
end

function detectoutbreak!(outbreakvec, incvec, avgvec, threshold)
    @. outbreakvec = ifelse(incvec >= threshold || avgvec >= threshold, 1, 0)

    return nothing
end

function calculate_test_positivity(
    true_positive_vec, total_positive_vec, alert_vec, agg_days
)
    @views outvec = zeros(Float64, length(true_positive_vec) ÷ agg_days, 2)
    @inbounds for i in axes(outvec, 1)
        start_ind = 1 + (i - 1) * agg_days
        end_ind = start_ind + (agg_days - 1)

        @views total_positive_sum = sum(total_positive_vec[start_ind:end_ind])
        @views true_positive_sum = sum(true_positive_vec[start_ind:end_ind])
        @views num_outbreak_days = sum(alert_vec[start_ind:end_ind])
        agg_outbreak_status = num_outbreak_days >= agg_days / 2 ? 1 : 0

        outvec[i, 1] = true_positive_sum / total_positive_sum
        outvec[i, 2] = agg_outbreak_status
    end
    return outvec
end

function calculate_OutbreakThresholdChars(testarr, infecarr, thresholds_vec)
    OT_chars = map(axes(infecarr, 3)) do sim
        dailychars = calculate_daily_detection_characteristics(
            @view(testarr[:, 7, sim]), @view(infecarr[:, 3, sim])
        )
        alertrle = rle(@view(testarr[:, 7, sim]))
        outbreakbounds = thresholds_vec[sim][
            (@view(thresholds_vec[sim][:, 4]) .== 1), 1:3
        ]
        alertbounds = calculate_outbreak_thresholds(alertrle; ncols = 2)

        detectionchars = calculate_outbreak_detection_characteristics(
            outbreakbounds, alertbounds
        )

        first_matched_bounds = filter_first_matched_bounds(
            detectionchars.matched_bounds
        )
        delay_vec = calculate_delay_vec(first_matched_bounds)

        cases_before_alert_vec, caseperc_before_alert,
        cases_after_alert_vec, caseperc_after_alert_vec = calculate_cases_before_after_alert(
            @view(infecarr[:, 1, sim]), first_matched_bounds, delay_vec
        )

        unavoidable_cases = calculate_unavoidable_cases(
            detectionchars.missed_outbreak_size, cases_before_alert_vec
        )

        avoidable_cases = calculate_avoidable_cases(cases_after_alert_vec)

        OutbreakThresholdChars(
            dailychars...,
            detectionchars...,
            delay_vec,
            cases_before_alert_vec,
            caseperc_before_alert,
            cases_after_alert_vec,
            caseperc_after_alert_vec,
            unavoidable_cases,
            avoidable_cases,
        )
    end

    return StructArray(OT_chars)
end

function calculate_unavoidable_cases(
    missed_outbreak_vec, cases_before_alert_vec
)
    return sum(missed_outbreak_vec) + sum(cases_before_alert_vec)
end

function calculate_avoidable_cases(cases_after_alert_vec)
    return sum(cases_after_alert_vec)
end

function calculate_daily_detection_characteristics(testvec, infecvec)
    crosstab = freqtable(testvec, infecvec)

    tp = freqtable_error_default_zero(crosstab, 1, 1)
    fp = freqtable_error_default_zero(crosstab, 1, 0)
    tn = freqtable_error_default_zero(crosstab, 0, 0)
    fn = freqtable_error_default_zero(crosstab, 0, 1)

    sens = tp / (tp + fn)
    spec = tn / (tn + fp)

    ppv = tp / (tp + fp)
    npv = tn / (tn + fn)

    return crosstab, tp, tn, fp, fn, sens, spec, ppv, npv
end

# TODO: write tests that check the function works when crosstab produces a 2x2 table,
# a 2x1, and a 1x1 table
function freqtable_error_default_zero() end

function freqtable_error_default_zero(
    freqtable, testing_class::T, actual_class::T
) where {T<:Integer}
    return try
        freqtable[Name(testing_class), Name(actual_class)]
    catch
        0
    end
end

function calculate_noutbreaks(outbreakrle)
    return length(findall(==(1), outbreakrle[1]))
end

function calculate_outbreak_detection_characteristics(
    outbreakbounds, alertbounds
)
    filtered_matched_bounds, periodssum_vec, alerts_per_outbreak_vec = match_outbreak_detection_bounds(
        outbreakbounds, alertbounds
    )

    noutbreaks = size(outbreakbounds, 1)
    nalerts = size(alertbounds, 1)

    detected_outbreak_size = periodssum_vec[(alerts_per_outbreak_vec .> 0)]
    missed_outbreak_size = periodssum_vec[(alerts_per_outbreak_vec .== 0)]

    n_true_outbreaks_detected = length(
        Set(@view(filtered_matched_bounds[:, 1]))
    )

    n_missed_outbreaks = noutbreaks - n_true_outbreaks_detected

    n_correct_alerts = size(filtered_matched_bounds, 1)

    n_false_alerts = sum(nalerts - n_correct_alerts)

    perc_true_outbreaks_detected = n_true_outbreaks_detected / noutbreaks
    perc_true_outbreaks_missed = n_missed_outbreaks / noutbreaks
    falsealert_trueoutbreak_prop = n_false_alerts / noutbreaks
    correctalert_trueoutbreak_prop = n_correct_alerts / noutbreaks # c.f. sensitivity

    trueoutbreak_alerts_prop = noutbreaks / nalerts
    outbreaksmissed_alerts_prop = n_missed_outbreaks / nalerts
    perc_alerts_false = n_false_alerts / nalerts
    perc_alerts_correct = n_correct_alerts / nalerts # c.f. PPV

    accuracy = NaNMath.mean([perc_true_outbreaks_detected, perc_alerts_correct])

    return (
        accuracy = accuracy,
        matched_bounds = filtered_matched_bounds,
        noutbreaks = noutbreaks,
        nalerts = nalerts,
        detected_outbreak_size = detected_outbreak_size,
        missed_outbreak_size = missed_outbreak_size,
        n_true_outbreaks_detected = n_true_outbreaks_detected,
        n_missed_outbreaks = n_missed_outbreaks,
        n_correct_alerts = n_correct_alerts,
        n_false_alerts = n_false_alerts,
        alertsperoutbreak = alerts_per_outbreak_vec,
        periodsumvec = periodssum_vec,
        perc_true_outbreaks_detected = perc_true_outbreaks_detected,
        perc_true_outbreaks_missed = perc_true_outbreaks_missed,
        falsealert_trueoutbreak_prop = falsealert_trueoutbreak_prop,
        correctalert_trueoutbreak_prop = correctalert_trueoutbreak_prop,
        trueoutbreak_alerts_prop = trueoutbreak_alerts_prop,
        outbreaksmissed_alerts_prop = outbreaksmissed_alerts_prop,
        perc_alerts_false = perc_alerts_false,
        perc_alerts_correct = perc_alerts_correct,
    )
end

function match_outbreak_detection_bounds(outbreakbounds, alertbounds)
    all_matched_bounds = zeros(
        Int64, size(outbreakbounds, 1) + size(alertbounds, 1), 5
    )
    alerts_per_outbreak_vec = zeros(Int64, size(outbreakbounds, 1))
    periodssum_vec = zeros(Int64, size(outbreakbounds, 1))

    outbreak_number = 1
    alert_rownumber = 1
    for (outbreak_number, (outbreaklower, outbreakupper, periodsum)) in
        pairs(eachrow(outbreakbounds))
        periodssum_vec[outbreak_number] = periodsum
        for (alertlower, alertupper) in
            eachrow(@view(alertbounds[alert_rownumber:end, :]))
            if alertlower > outbreakupper
                break
            end
            if alertlower >= outbreaklower
                all_matched_bounds[alert_rownumber, :] .= outbreaklower,
                outbreakupper, alertlower, alertupper, periodsum

                alerts_per_outbreak_vec[outbreak_number] += 1

                alert_rownumber += 1
                continue
            end
            if alertlower <= outbreaklower &&
                alertupper > outbreaklower
                all_matched_bounds[alert_rownumber, :] .= outbreaklower,
                outbreakupper, alertlower, alertupper, periodsum

                alerts_per_outbreak_vec[outbreak_number] += 1

                alert_rownumber += 1
                continue
            end
        end
        outbreak_number += 1
    end

    filtered_matched_bounds = @view(
        all_matched_bounds[(all_matched_bounds[:, 2] .> 0), :]
    )
    return filtered_matched_bounds, periodssum_vec, alerts_per_outbreak_vec
end

function calculate_delay_vec(first_matchedbounds)
    return @views first_matchedbounds[:, 3] .- first_matchedbounds[:, 1]
end

function filter_first_matched_bounds(matchedbounds)
    indices = calculate_first_matched_bounds_index(matchedbounds)
    return matchedbounds[indices, :]
end

function calculate_first_matched_bounds_index(matchedbounds)
    return map(
        outbreaklower -> findfirst(
            isequal(outbreaklower),
            @view(matchedbounds[:, 1])
        ),
        unique(@view(matchedbounds[:, 1])),
    )
end

function calculate_cases_before_after_alert(
    incvec, first_matchedbounds, delay_vec
)
    casesaftervec = zeros(Int64, length(delay_vec))
    casesafterpercvec = zeros(Float64, length(casesaftervec))
    casesbeforevec = zeros(Int64, length(delay_vec))
    casesbeforepercvec = zeros(Float64, length(casesaftervec))
    calculate_cases_before_after_alert!(
        casesbeforevec,
        casesbeforepercvec,
        casesaftervec,
        casesafterpercvec,
        incvec,
        first_matchedbounds,
        delay_vec,
    )
    return casesbeforevec, casesbeforepercvec, casesaftervec, casesafterpercvec
end

function calculate_cases_before_after_alert!(
    casesbeforevec,
    casesbeforepercvec,
    casesaftervec,
    casesafterpercvec,
    invec,
    first_matchedbounds,
    delay_vec,
)
    for (
        alertnumber,
        (outbreaklower, outbreakupper, alertlower, alertupper, periodsum),
    ) in
        pairs(eachrow(first_matchedbounds))
        if delay_vec[alertnumber] <= 0
            casesaftervec[alertnumber, 1] = periodsum
            casesbeforevec[alertnumber, 1] = 0
        else
            @views casesaftervec[alertnumber] = sum(
                invec[alertlower:outbreakupper]
            )
            casesbeforevec[alertnumber] = periodsum - casesaftervec[alertnumber]
        end
        casesbeforepercvec[alertnumber] =
            casesbeforevec[alertnumber] / periodsum
        casesafterpercvec[alertnumber] = casesaftervec[alertnumber] / periodsum
    end
    return nothing
end

# end
