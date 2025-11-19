export objective_function, calculate_ensemble_objective_metric,
    calculate_outbreak_detection_accuracy

"""
    objective_function(alert_threshold_vec, inputs)

Objective function for threshold optimization.

Returns 1 - mean_accuracy (for minimization).
"""
function objective_function(
        alert_threshold_vec,
        inputs,
    )
    @assert length(alert_threshold_vec) == 1

    UnPack.@unpack ensemble_inc_arr,
        noise_array,
        outbreak_detection_specification,
        individual_test_specification,
        thresholds_vec,
        accuracy_function = inputs

    outbreak_detection_specification = OutbreakDetectionSpecification(
        alert_threshold_vec[1],
        outbreak_detection_specification.moving_average_lag,
        outbreak_detection_specification.percent_visit_clinic,
        outbreak_detection_specification.percent_clinic_tested,
        outbreak_detection_specification.alert_method.method_name,
    )

    testarr = create_testing_arrs(
        ensemble_inc_arr,
        noise_array,
        outbreak_detection_specification,
        individual_test_specification,
    )[1]

    objective = calculate_ensemble_objective_metric(
        accuracy_function, testarr, ensemble_inc_arr, thresholds_vec
    )

    return objective
end

"""
    calculate_ensemble_objective_metric(accuracy_function, testarr, 
                                        infecarr, thresholds_vec)

Calculate ensemble-level objective metric for optimization.
"""
function calculate_ensemble_objective_metric(
        accuracy_function, testarr, infecarr, thresholds_vec
    )
    mean_accuracy = map(axes(infecarr, 3)) do sim
        dailychars = calculate_daily_detection_characteristics(
            @view(testarr[:, 6, sim]), @view(infecarr[:, 3, sim])
        )
        alertrle = StatsBase.rle(@view(testarr[:, 6, sim]))
        outbreakbounds = thresholds_vec[sim]
        alertbounds = calculate_outbreak_thresholds(alertrle; ncols = 3)

        accuracy = calculate_outbreak_detection_accuracy(
            accuracy_function, outbreakbounds, alertbounds
        )
    end

    return 1 - NaNMath.mean(mean_accuracy)
end

"""
    calculate_outbreak_detection_accuracy(accuracy_function, 
                                          outbreakbounds, alertbounds)

Calculate outbreak detection accuracy for a single simulation.
"""
function calculate_outbreak_detection_accuracy(
        accuracy_function, outbreakbounds, alertbounds
    )
    filtered_matched_bounds = match_outbreak_detection_bounds(
        outbreakbounds, alertbounds
    )[1]

    noutbreaks = size(outbreakbounds, 1)
    nalerts = size(alertbounds, 1)

    n_true_outbreaks_detected = length(
        Set(@view(filtered_matched_bounds[:, 1]))
    )
    n_correct_alerts = size(filtered_matched_bounds, 1)

    perc_true_outbreaks_detected = n_true_outbreaks_detected / noutbreaks
    perc_alerts_correct = n_correct_alerts / nalerts # c.f. PPV

    accuracy = accuracy_function(
        perc_alerts_correct, perc_true_outbreaks_detected
    )

    return accuracy
end
