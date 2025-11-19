using StatsBase: StatsBase
using FreqTables: FreqTables
using StructArrays: StructArrays
using NaNMath: NaNMath

export calculate_OutbreakThresholdChars, calculate_daily_detection_characteristics,
    calculate_outbreak_detection_characteristics, calculate_noutbreaks,
    calculate_n_outbreak_tests, filter_first_matched_bounds,
    calculate_first_matched_bounds_index, calculate_cases_before_after_alert!,
    calculate_cases_before_after_alert, calculate_f_beta_score, arithmetic_mean

"""
    calculate_OutbreakThresholdChars(testarr, infecarr, thresholds_vec, noise_means)

Calculate comprehensive outbreak detection characteristics.
"""
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

"""
    calculate_n_outbreak_tests(infectious_tested_vec, noise_tested_vec, outbreakbounds)

Calculate total number of tests during outbreak periods.
"""
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

"""
    calculate_daily_detection_characteristics(testvec, infecvec)

Calculate daily-level detection characteristics (sensitivity, specificity, PPV, NPV).
"""
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
    ) where {T <: Integer}
    return try
        freqtable[FreqTables.Name(testing_class), FreqTables.Name(actual_class)]
    catch
        0
    end
end

"""
    calculate_noutbreaks(outbreakrle)

Calculate number of outbreaks from run-length encoding.
"""
function calculate_noutbreaks(outbreakrle)
    return length(findall(==(1), outbreakrle[1]))
end

"""
    calculate_outbreak_detection_characteristics(outbreakbounds, alertbounds)

Calculate outbreak-level detection characteristics.
"""
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

Generic formula to calculate the arithmetic mean. Used for the `accuracy` 
measure for outbreak detection.

  - Precision = PPV (% alerts that are correct)
  - Recall = sensitivity (% of outbreaks detected)

Implement as `arithmetic_mean(perc_alerts_correct, perc_true_outbreaks_detected)`
"""
function arithmetic_mean(precision, recall)
    return NaNMath.mean([precision, recall])
end

"""
    calculate_f_beta_score(precision, recall; beta = 1)

Generic formula to calculate the F-score. When beta=1, calculates the F1 score, 
which weights precision and recall equally (the harmonic mean). beta=2 weights 
precision more than recall, and beta=0.5 weights recall more.

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

"""
    filter_first_matched_bounds(matchedbounds)

Filter to keep only the first alert for each outbreak.
"""
function filter_first_matched_bounds(matchedbounds)
    indices = calculate_first_matched_bounds_index(matchedbounds)
    return matchedbounds[indices, :]
end

"""
    calculate_first_matched_bounds_index(matchedbounds)

Calculate indices of first matched bounds for each outbreak.
"""
function calculate_first_matched_bounds_index(matchedbounds)
    return map(
        outbreaklower -> findfirst(
            isequal(outbreaklower),
            @view(matchedbounds[:, 1])
        ),
        unique(@view(matchedbounds[:, 1])),
    )
end

"""
    calculate_cases_before_after_alert(incvec, first_matchedbounds, delay_vec)

Calculate cases before and after alert for each outbreak.
"""
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

"""
    calculate_cases_before_after_alert!(casesbeforevec, casesbeforepercvec, 
                                        casesaftervec, casesafterpercvec, 
                                        invec, first_matchedbounds, delay_vec)

In-place calculation of cases before and after alert.
"""
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
