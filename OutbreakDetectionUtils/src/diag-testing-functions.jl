using StatsBase: StatsBase
using FreqTables: FreqTables
using StructArrays: StructArrays
using NaNMath: NaNMath
using Match: Match

function create_testing_arrs(
    incarr,
    noisearr,
    outbreak_detect_spec::OutbreakDetectionSpecification,
    individual_test_spec::IndividualTestSpecification,
)
    testarr = zeros(Int64, size(incarr, 1), 7, size(incarr, 3))
    test_movingavg_arr = zeros(Int64, size(incarr, 1), size(incarr, 3))

    create_testing_arrs!(
        testarr,
        test_movingavg_arr,
        incarr,
        noisearr,
        outbreak_detect_spec.alert_method.method_name,
        outbreak_detect_spec.alert_threshold,
        outbreak_detect_spec.moving_average_lag,
        outbreak_detect_spec.percent_tested,
        individual_test_spec.test_result_lag,
        individual_test_spec.sensitivity,
        individual_test_spec.specificity,
    )

    return testarr, test_movingavg_arr
end

function create_testing_arrs!(
    testarr,
    test_movingavg_arr,
    incarr,
    noisearr,
    alert_method,
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
            @view(testarr[:, 2, sim]), @view(noisearr[:, sim]), perc_tested
        )

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
            @view(test_movingavg_arr[:, sim]),
            @view(testarr[:, 5, sim]),
            moveavglag,
        )

        detectoutbreak_args = Match.@match alert_method begin
            "movingavg" => (
                @view(testarr[:, 6, sim]),
                @view(test_movingavg_arr[:, sim]),
                alertthreshold,
            )
            "dailythreshold_movingavg" => (
                @view(testarr[:, 6, sim]),
                @view(testarr[:, 5, sim]),
                @view(test_movingavg_arr[:, sim]),
                alertthreshold,
            )
        end

        # TOTAL Test positive individuals trigger outbreak response
        detectoutbreak!(detectoutbreak_args...)

        # Triggered outbreak equal to actual outbreak status
        @. testarr[:, 7, sim] =
            @view(testarr[:, 6, sim]) == @view(incarr[:, 3, sim])
    end

    return nothing
end

function calculate_tested!(outvec, invec, perc_tested)
    @. outvec = round(invec * perc_tested)
end

function calculate_positives(
    type_positive_function!, tested_vec, tlength, lag, testcharacteristic
)
    outvec = zeros(Int64, tlength)
    type_positive_function!(
        outvec, tested_vec, tlength, lag, testcharacteristic
    )
    return outvec
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
    return StatsBase.mean(@view(invec[moveavg_daystart:day]))
end

function calculate_int_daily_movingavg(invec, day, avglag)
    @inline moveavg_daystart = calculate_daily_movingavg_startday(day, avglag)
    return Int64(round(StatsBase.mean(@view(invec[moveavg_daystart:day]))))
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

"""
    detectoutbreak!(outbreakvec, incvec, avgvec, threshold)

Determines whether an outbreak has been detected based on an infection daily and moving average timeseries and a threshold.
`incvec` should be the daily incidence timeseries, and `avgvec` should be the moving average of the daily incidence timeseries.
"""
function detectoutbreak!(outbreakvec, incvec, avgvec, threshold)
    @. outbreakvec = ifelse(incvec >= threshold || avgvec >= threshold, 1, 0)

    return nothing
end

function detectoutbreak(infectionsvec, threshold)
    outbreak = zeros(Int64, length(infectionsvec))

    detectoutbreak!(outbreak, infectionsvec, threshold)

    return outbreak
end

"""
    detectoutbreak!(outbreakvec, infectionsvec, threshold)

Determines whether an outbreak has been detected based on an infection timeseries and a threshold.
The `infectionsvec` can either be a vector of daily infections or a vector of the moving average of daily infections.
"""
function detectoutbreak!(outbreakvec, infectionsvec, threshold)
    @. outbreakvec = ifelse(infectionsvec >= threshold, 1, 0)

    return nothing
end

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

function calculate_OutbreakThresholdChars(
    testarr, infecarr, thresholds_vec, noise_means
)
    OT_chars = map(axes(infecarr, 3)) do sim
        dailychars = calculate_daily_detection_characteristics(
            @view(testarr[:, 6, sim]), @view(infecarr[:, 3, sim])
        )
        alertrle = StatsBase.rle(@view(testarr[:, 6, sim]))
        outbreakbounds = thresholds_vec[sim]
        alertbounds = calculate_outbreak_thresholds(alertrle; ncols = 3)
        calculate_outbreak_duration!(alertbounds)

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

        n_outbreak_cases = unavoidable_cases + avoidable_cases

        n_tests = calculate_n_tests(
            @view(testarr[:, 1, sim]), @view(testarr[:, 2, sim])
        )

        n_outbreak_tests = calculate_n_outbreak_tests(
            @view(testarr[:, 1, sim]), @view(testarr[:, 2, sim]),
            outbreakbounds,
        )

        mean_noise_incidence_ratio =
            noise_means.mean_noise / StatsBase.mean(@view(infecarr[:, 1, sim]))

        proportion_timeseries_in_outbreak = calculate_proportion_timeseries_in_outbreak(
            @view(infecarr[:, 3, sim])
        )

        proportion_timeseries_in_alert = calculate_proportion_timeseries_in_outbreak(
            @view(testarr[:, 6, sim])
        )

        alert_outbreak_timeseries_prop_diff =
            proportion_timeseries_in_alert - proportion_timeseries_in_outbreak

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
            n_outbreak_cases,
            n_tests,
            n_outbreak_tests,
            mean_noise_incidence_ratio,
            noise_means.mean_poisson_noise,
            noise_means.poisson_noise_prop,
            proportion_timeseries_in_outbreak,
            proportion_timeseries_in_alert,
            alert_outbreak_timeseries_prop_diff,
        )
    end

    return StructArrays.StructArray(OT_chars)
end

function calculate_proportion_timeseries_in_outbreak(outbreak_status_vec)
    return sum(outbreak_status_vec) / length(outbreak_status_vec)
end

function calculate_proportion_timeseries_in_outbreak(
    outbreak_bounds_periodsum_vec, time_parameters
)
    return sum(outbreak_bounds_periodsum_vec) / time_parameters.tlength
end

function calculate_n_tests(infectious_tested_vec, noise_tested_vec)
    return sum(infectious_tested_vec) + sum(noise_tested_vec)
end

function calculate_n_outbreak_tests(
    infectious_tested_vec, noise_tested_vec, outbreakbounds
)
    tested = 0
    for (lower, upper) in eachrow(@view(outbreakbounds[:, 1:2]))
        @views tested += sum(infectious_tested_vec[lower:upper])
        @views tested += sum(noise_tested_vec[lower:upper])
    end
    return tested
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
    crosstab = FreqTables.freqtable(testvec, infecvec)

    tp = freqtable_error_default_zero(crosstab, 1, 1)
    fp = freqtable_error_default_zero(crosstab, 1, 0)
    tn = freqtable_error_default_zero(crosstab, 0, 0)
    fn = freqtable_error_default_zero(crosstab, 0, 1)

    sens = tp / (tp + fn)
    spec = tn / (tn + fp)

    ppv = tp / (tp + fp)
    npv = tn / (tn + fn)

    return sens, spec, ppv, npv
end

# TODO: write tests that check the function works when crosstab produces a 2x2 table,
# a 2x1, and a 1x1 table
function freqtable_error_default_zero() end

function freqtable_error_default_zero(
    freqtable, testing_class::T, actual_class::T
) where {T<:Integer}
    return try
        freqtable[FreqTables.Name(testing_class), FreqTables.Name(actual_class)]
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
    filtered_matched_bounds, outbreak_dur_vec, alert_dur_vec, periodssum_vec, alerts_per_outbreak_vec = match_outbreak_detection_bounds(
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

    accuracy = arithmetic_mean(
        perc_alerts_correct, perc_true_outbreaks_detected
    )
    f1_score = calculate_f_beta_score(
        perc_alerts_correct, perc_true_outbreaks_detected
    )

    return (
        accuracy = accuracy,
        f1_score = f1_score,
        matched_bounds = filtered_matched_bounds,
        noutbreaks = noutbreaks,
        nalerts = nalerts,
        outbreak_duration_vec = outbreak_dur_vec,
        alert_duration_vec = alert_dur_vec,
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
    periodssum_vec = outbreakbounds[:, 4]
    outbreak_dur_vec = outbreakbounds[:, 3]
    alert_dur_vec = alertbounds[:, 3]

    outbreak_number = 1
    alert_rownumber = 1
    for (outbreak_number, (outbreaklower, outbreakupper, periodsum)) in
        pairs(eachrow(@view(outbreakbounds[:, [1, 2, 4]])))
        for (alertlower, alertupper) in
            eachrow(@view(alertbounds[alert_rownumber:end, 1:2]))
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
    return filtered_matched_bounds,
    outbreak_dur_vec, alert_dur_vec, periodssum_vec,
    alerts_per_outbreak_vec
end

"""
    arithmetic_mean(precision, recall)

Generic formula to calculate the arithmetic mean. Used for the `accuracy` measure for outbreak detection.

  - Precision = PPV (% alerts that are correct)
  - Recall = sensitivity (% of outbreaks detected)

Implement as `arithmetic_mean(perc_alerts_correct, perc_true_outbreaks_detected)`
"""
function arithmetic_mean(precision, recall)
    return NaNMath.mean([precision, recall])
end

"""
    calculate_f_beta_score(precision, recall; beta = 1)

Generic formula to calculate the F-score. When beta=1, calculates the F1 score, which weights precision and recall equally (the harmonic mean). beta=2 weights precision more than recall, and beta=0.5 weights recall more.

  - Precision = PPV (% alerts that are correct)
  - Recall = sensitivity (% of outbreaks detected)

Implement as `calculate_f_beta_score(perc_alerts_correct, perc_true_outbreaks_detected; beta = 1)`
"""
function calculate_f_beta_score(
    precision, recall; beta = 1
)
    f_beta =
        (1 + beta^2) * (precision * recall) / ((beta^2 * precision) + recall)

    if isnan(f_beta)
        f_beta = 0.0
    end

    return f_beta
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
