function classify_outbreak_detection(
        optimization_scenario::OptimizationScenario,
        emergent_ews_metrics::StructVector{EWSMetrics},
        null_ews_metrics::StructVector{EWSMetrics},
    )
    @unpack ews_metric,
        ews_threshold_window,
        ensemble_specification,
        threshold_quantile,
        consecutive_thresholds = optimization_scenario
    @unpack burnin = ensemble_specification.time_parameters

    ews_metric_symbol = Symbol(ews_metric)

    nsims = length(emergent_ews_metrics)

    true_positives = 0

    for sim in eachindex(emergent_ews_metrics)
        # Use pre-computed EWS metrics
        ews_vals = emergent_ews_metrics[sim]
        null_ews_vals = null_ews_metrics[sim]

        # Check threshold exceedances
        exceeds_threshold = exceeds_ews_threshold(
            ews_vals,
            ews_metric_symbol,
            ews_threshold_window,
            threshold_quantile,
            burnin,
        )

        detection_index = calculate_ews_trigger_index(
            exceeds_threshold,
            consecutive_thresholds,
        )

        null_exceeds_threshold = exceeds_ews_threshold(
            null_ews_vals,
            ews_metric_symbol,
            ews_threshold_window,
            threshold_quantile,
            burnin,
        )

        null_detection_index = calculate_ews_trigger_index(
            null_exceeds_threshold,
            consecutive_thresholds,
        )

        # Update counts
        if Try.isok(detection_index)
            true_positives += 1
        end
        if Try.iserr(null_detection_index)
            true_negatives += 1
        end
    end

    return AlertClassificationResults(
        true_positives,
        n_null_sims - true_negatives,
        n_emergent_sims - true_positives,
        n_sims,
    )
end
