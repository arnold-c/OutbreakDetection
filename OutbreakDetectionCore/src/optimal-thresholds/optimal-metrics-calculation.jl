export calculate_optimal_results

"""
    calculate_optimal_results(
        alert_threshold,
        accuracy_metric,
        test_positives_container,
        ensemble_simulation,
        outbreak_thresholds
    )

Calculate comprehensive detection performance metrics at a given alert threshold.

This function evaluates outbreak detection performance across all simulations in an
ensemble by generating alerts at the specified threshold and comparing them against
known outbreak periods. It computes a wide range of metrics including accuracy,
detection delays, unavoidable cases, and temporal characteristics.

# Arguments

  - `alert_threshold::Float64`: Threshold value for generating outbreak alerts
  - `accuracy_metric::AccuracyMetric`: Metric for calculating accuracy (e.g., F1Score, PPV, Sensitivity)
  - `test_positives_container::StructVector{TestPositiveContainer}`: Container with test positive data for each simulation, structured according to the alert method
  - `ensemble_simulation::StructVector{SEIRRun}`: Collection of SEIR simulation runs containing incidence trajectories
  - `outbreak_thresholds::StructVector{OutbreakThresholds}`: Pre-computed outbreak threshold bounds for each simulation

# Returns

A named tuple containing vectors of metrics for each simulation:

  - `accuracies`: Detection accuracy based on the specified metric
  - `proportion_alerts_correct`: Positive predictive value (PPV) - proportion of alerts that correctly identified outbreaks
  - `proportion_outbreaks_detected`: Sensitivity - proportion of outbreaks that were detected
  - `detection_delays`: Mean delay (in time steps) between outbreak start and alert generation
  - `unavoidable_cases`: Number of cases that occurred before detection (unavoidable given detection delay)
  - `alert_durations`: Mean duration of alert periods
  - `outbreak_durations`: Mean duration of outbreak periods
  - `proportion_timeseries_in_alert`: Proportion of the time series spent in alert state
  - `proportion_timeseries_in_outbreak`: Proportion of the time series spent in outbreak state

# Details

For each simulation, the function:

  1. Generates binary alert vector based on the threshold and alert method
  2. Identifies alert periods using run-length encoding
  3. Matches alert periods to outbreak periods
  4. Calculates PPV (proportion of alerts that are correct)
  5. Calculates sensitivity (proportion of outbreaks detected)
  6. Computes accuracy using the specified metric
  7. Calculates detection delays and unavoidable cases
  8. Computes temporal characteristics (durations and proportions)
"""
function calculate_optimal_results(
        alert_threshold::Float64,
        accuracy_metric::AccuracyMetric,
        test_positives_container::StructVector{TestPositiveContainer},
        ensemble_simulation::StructVector{SEIRRun},
        outbreak_thresholds::StructVector{OutbreakThresholds},
    )
    nsims = length(outbreak_thresholds)

    # Accumulators for metrics
    accuracies = Vector{Float64}(undef, nsims)
    proportion_alerts_correct = similar(accuracies)
    proportion_outbreaks_detected = similar(accuracies)
    detection_delays = Vector{Vector{Int64}}(undef, nsims)
    unavoidable_cases = similar(detection_delays)
    alert_durations = similar(detection_delays)
    outbreak_durations = similar(alert_durations)
    proportion_timeseries_in_alert = similar(accuracies)
    proportion_timeseries_in_outbreak = similar(accuracies)

    for sim in eachindex(accuracies)
        # Generate alerts based on threshold and alert method
        alert_vec = generate_alerts(
            test_positives_container[sim],
            alert_threshold
        )

        # Get pre-computed outbreak bounds for this simulation
        outbreak_bounds = outbreak_thresholds[sim]

        # Get alert bounds using RLE
        alert_rle = StatsBase.rle(alert_vec)
        alert_bounds = calculate_above_threshold_bounds(alert_rle)

        matched_outbreak_threshold_indices = match_outbreak_detection_bounds(
            outbreak_bounds,
            alert_bounds
        )

        ppv = calculate_ppv(matched_outbreak_threshold_indices)
        proportion_alerts_correct[sim] = ppv

        sensitivity = calculate_sensitivity(matched_outbreak_threshold_indices)
        proportion_outbreaks_detected[sim] = sensitivity

        accuracies[sim] = calculate_accuracy(
            accuracy_metric,
            ppv,
            sensitivity
        )

        detection_delays[sim] = calculate_detection_delay(
            matched_outbreak_threshold_indices,
            outbreak_bounds,
            alert_bounds,
        )

        unavoidable_cases[sim] = calculate_unavoidable_cases(
            matched_outbreak_threshold_indices,
            outbreak_bounds,
            alert_bounds,
            ensemble_simulation[sim].incidence
        )

        alert_durations[sim] = get_alert_duration(alert_bounds)
        outbreak_durations[sim] = get_outbreak_duration(outbreak_bounds)

        proportion_timeseries_in_alert[sim] = calculate_proportion_timeseries_in_alert(
            alert_bounds,
            length(ensemble_simulation[sim].incidence)
        )
        proportion_timeseries_in_outbreak[sim] = calculate_proportion_timeseries_in_outbreak(
            outbreak_bounds,
            length(ensemble_simulation[sim].incidence)
        )
    end

    return (;
        accuracies,
        proportion_alerts_correct,
        proportion_outbreaks_detected,
        detection_delays,
        unavoidable_cases,
        alert_durations,
        outbreak_durations,
        proportion_timeseries_in_alert,
        proportion_timeseries_in_outbreak,
    )
end
