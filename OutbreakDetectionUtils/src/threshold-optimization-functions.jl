using NaNMath: NaNMath

function calculate_ensemble_objective_metric(
    testarr, infecarr, thresholds_vec, noise_means
)
    mean_accuracy = mapreduce(NaNMath.mean, axes(infecarr, 3)) do sim
        dailychars = calculate_daily_detection_characteristics(
            @view(testarr[:, 6, sim]), @view(infecarr[:, 3, sim])
        )
        alertrle = StatsBase.rle(@view(testarr[:, 6, sim]))
        outbreakbounds = thresholds_vec[sim]
        alertbounds = calculate_outbreak_thresholds(alertrle; ncols = 3)
        # calculate_outbreak_duration!(alertbounds)

        accuracy = calculate_outbreak_detection_accuracy(
            outbreakbounds, alertbounds
        )
    end

	return 1 - mean_accuracy
end

function calculate_outbreak_detection_accuracy(
    outbreakbounds, alertbounds
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

    accuracy = NaNMath.mean([perc_true_outbreaks_detected, perc_alerts_correct])

    return accuracy
end
