function plot_all_threshold_comparisons(percent_clinic_tested, base_parameters)
    @unpack sensitivity_vec,
    specificity_vec,
    detectthreshold_vec,
    ensemble_specification,
    noise_specification,
    outbreak_specification,
    moving_avg_detection_lag,
    test_result_lag,
    percent_visit_clinic = base_parameters

    ensemble_scenario_spec_vec = Vector{ScenarioSpecification}(
        undef,
        length(sensitivity_vec) * length(detectthreshold_vec) +
        length(detectthreshold_vec),
    )
    ensemble_chars_vec = Vector(
        undef, length(ensemble_scenario_spec_vec)
    )

    outbreak_detect_spec_vec = map(
        threshold -> OutbreakDetectionSpecification(
            threshold,
            moving_avg_detection_lag,
            percent_visit_clinic,
            percent_clinic_tested,
            test_result_lag,
        ),
        detectthreshold_vec,
    )
    clinical_case_outbreak_detect_spec_vec = map(
        threshold -> OutbreakDetectionSpecification(
            threshold,
            moving_avg_detection_lag,
            percent_visit_clinic,
            1.0,
            test_result_lag,
        ),
        detectthreshold_vec,
    )

    for (i, ((sens, spec), outbreak_detect_spec)) in enumerate(
        Iterators.product(
            zip(sensitivity_vec, specificity_vec),
            outbreak_detect_spec_vec
        ),
    )
        ind_test_spec = IndividualTestSpecification(sens, spec)

        ensemble_scenario_spec = ScenarioSpecification(
            ensemble_specification,
            outbreak_specification,
            noise_specification,
            outbreak_detect_spec,
            ind_test_spec,
        )

        ensemble_scenario_spec_vec[i] = ensemble_scenario_spec
    end

    ensemble_scenario_spec_vec[(end - length(clinical_case_outbreak_detect_spec_vec) + 1):end] .= create_combinations_vec(
        ScenarioSpecification,
        (
            [ensemble_specification],
            [outbreak_specification],
            [noise_specification],
            clinical_case_outbreak_detect_spec_vec,
            # TODO: update this to calculate for all detection thresholds
            [IndividualTestSpecification(1.0, 0.0)],
        ),
    )

    prog = Progress(length(ensemble_scenario_spec_vec))
    @floop for (i, ensemble_scenario_spec) in pairs(ensemble_scenario_spec_vec)
        ensemble_chars_file = get_ensemble_file(ensemble_scenario_spec)

        ensemble_chars_vec[i] = (
            OT_chars = ensemble_chars_file["OT_chars"],
            outbreak_detect_spec = ensemble_scenario_spec.outbreak_detection_specification,
            ind_test_spec = ensemble_scenario_spec.individual_test_specification,
        )
        next!(prog)
    end

    sort!(
        ensemble_chars_vec;
        by = x -> (
            x.outbreak_detect_spec.detection_threshold,
            x.ind_test_spec.specificity,
        ),
    )

    compare_outbreak_sens_spec_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
            (
                char = :daily_sensitivity,
                label = "Sensitivity",
                color = (DAILY_SENSITIVITY_COLOR, 0.7),
            ),
            (
                char = :daily_specificity,
                label = "Specificity",
                color = (DAILY_SPECIFICITY_COLOR, 0.7),
            ),
        ];
        columnfacetchar_label = "Detection Threshold",
        bins = 0.0:0.01:1.01,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_sens_spec_plot.png"
        ),
        compare_outbreak_sens_spec_plot;
        resolution = (2200, 1200),
    )

    compare_outbreak_ppv_npv_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
            (
                char = :daily_ppv,
                label = "PPV",
                color = (DAILY_PPV_COLOR, 0.7),
            ),
            (
                char = :daily_npv,
                label = "NPV",
                color = (DAILY_NPV_COLOR, 0.7),
            ),
        ];
        columnfacetchar_label = "Detection Threshold",
        bins = 0.0:0.01:1.01,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_ppv_npv_plot.png"
        ),
        compare_outbreak_ppv_npv_plot;
        resolution = (2200, 1200),
    )

    compare_outbreak_detection_delays_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
        (
            char = :detectiondelays,
            color = (DETECTION_DELAY_COLOR, 1.0),
        )
    ];
        xlabel = "Detection Delay (days)",
        columnfacetchar_label = "Detection Threshold",
        binwidth = 5.0,
        meanlines = true,
        legend = false,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_detection_delays_plot.png",
        ),
        compare_outbreak_detection_delays_plot;
        resolution = (2200, 1200),
    )

    compare_outbreak_alert_per_outbreak_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
        (
            char = :n_alerts_per_outbreak,
            color = (N_ALERTS_PER_OUTBREAK_COLOR, 1.0),
        )
    ];
        xlabel = "Alerts per Outbreak",
        columnfacetchar_label = "Detection Threshold",
        binwidth = 1.0,
        meanlines = true,
        legend = false,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_alert_per_outbreak_plot.png",
        ),
        compare_outbreak_alert_per_outbreak_plot;
        resolution = (2200, 1200),
    )

    compare_outbreak_false_alerts_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
        (
            char = :n_false_alerts,
            color = (N_FALSE_ALERTS_COLOR, 1.0),
        )
    ];
        xlabel = "# False Alerts",
        columnfacetchar_label = "Detection Threshold",
        binwidth = 1.0,
        meanlines = true,
        legend = false,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_false_alerts_plot.png"
        ),
        compare_outbreak_false_alerts_plot;
        resolution = (2200, 1200),
    )

    compare_outbreak_number_alerts_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
        (
            char = :ndetectoutbreaks,
            color = (N_ALERTS_COLOR, 1.0),
        )
    ];
        xlabel = "# Alerts",
        columnfacetchar_label = "Detection Threshold",
        binwidth = 10.0,
        meanlines = true,
        legend = false,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_number_alerts_plot.png",
        ),
        compare_outbreak_number_alerts_plot;
        resolution = (2200, 1200),
    )

    compare_outbreak_numbers_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
        (
            char = :noutbreaks,
            color = (N_OUTBREAKS_COLOR, 1.0),
        )
    ];
        xlabel = "# Outbreaks",
        columnfacetchar_label = "Detection Threshold",
        binwidth = 1.0,
        meanlines = true,
        legend = false,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_numbers_plot.png"
        ),
        compare_outbreak_numbers_plot;
        resolution = (2200, 1200),
    )

    compare_outbreak_missed_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
        (
            char = :n_missed_outbreaks,
            color = (N_MISSED_OUTBREAKS_COLOR, 1.0),
        )
    ];
        xlabel = "# Missed Outbreaks",
        columnfacetchar_label = "Detection Threshold",
        binwidth = 1.0,
        meanlines = true,
        legend = false,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_missed_plot.png"
        ),
        compare_outbreak_missed_plot;
        resolution = (2200, 1200),
    )

    compare_outbreak_detect_missed_size_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
            (
                char = :detected_outbreak_size,
                label = "Size of Outbreaks Detected",
                color = (PERC_OUTBREAKS_DETECTED_COLOR, 0.7),
            ),
            (
                char = :missed_outbreak_size,
                label = "Size of Outbreaks Missed",
                color = (PERC_OUTBREAKS_MISSED_COLOR, 0.7),
            ),
        ];
        columnfacetchar_label = "Detection Threshold",
        binwidth = 200,
        normalization = :pdf,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_detect_missed_size_plot.png",
        ),
        compare_outbreak_detect_missed_size_plot;
        resolution = (2200, 1200),
    )

    compare_outbreak_true_outbreak_perc_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
            (
                char = :perc_true_outbreaks_detected,
                label = "Percent Outbreaks Detected",
                color = (PERC_OUTBREAKS_DETECTED_COLOR, 0.7),
            ),
            (
                char = :perc_true_outbreaks_missed,
                label = "Percent Outbreaks Missed",
                color = (PERC_OUTBREAKS_MISSED_COLOR, 0.7)),
        ];
        columnfacetchar_label = "Detection Threshold",
        binwidth = 0.02,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_true_outbreak_perc_plot.png",
        ),
        compare_outbreak_true_outbreak_perc_plot;
        resolution = (2200, 1200),
    )

    for i in eachindex(ensemble_chars_vec)
        perc_alerts_sum = sum(
            ensemble_chars_vec[i].OT_chars.perc_true_outbreaks_detected .+
            ensemble_chars_vec[i].OT_chars.perc_true_outbreaks_missed .- 1.0,
        )
        if perc_alerts_sum != 0.0
            @error "Warning. Sum of perc_alerts_false + perc_alerts_correct != 1.0, for i = $i. perc_alerts_sum = $perc_alerts_sum"
        end
    end

    compare_outbreak_alerts_perc_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
            (
                char = :perc_alerts_correct,
                label = "Percent Alerts\nThat Are Correct",
                color = (PERC_ALERTS_CORRECT_COLOR, 0.7),
            ),
            (
                char = :perc_alerts_false,
                label = "Percent Alerts\nThat Are False",
                color = (PERC_ALERTS_FALSE_COLOR, 0.7)),
        ];
        columnfacetchar_label = "Detection Threshold",
        bins = -0.01:0.02:1.01,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_alerts_perc_plot.png"
        ),
        compare_outbreak_alerts_perc_plot;
        resolution = (2200, 1200),
    )

    for i in eachindex(ensemble_chars_vec)
        perc_alerts_sum = sum(
            ensemble_chars_vec[i].OT_chars.perc_alerts_false .+
            ensemble_chars_vec[i].OT_chars.perc_alerts_correct .- 1.0,
        )
        if perc_alerts_sum != 0.0
            @error "Warning. Sum of perc_alerts_false + perc_alerts_correct != 1.0, for i = $i. perc_alerts_sum = $perc_alerts_sum,\nTimes no outbreaks detected = $(length(findall(==(0), ensemble_chars_vec[i].OT_chars.ndetectoutbreaks)))"
            nan_perc_alerts_sum = NaNMath.sum(
                ensemble_chars_vec[i].OT_chars.perc_alerts_false .+
                ensemble_chars_vec[i].OT_chars.perc_alerts_correct .- 1.0,
            )
            if nan_perc_alerts_sum != 0.0
                @error "Ignoring NaN values in the percentage of alerts doesn't correct the issue"
                continue
            end
            @info "Ignoring NaN values in the percentage of alerts does correct the issue"
        end
    end

    compare_outbreak_true_outbreak_alerts_perc_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
            (
                char = :perc_alerts_correct,
                label = "Percent Alerts\nThat Are Correct",
                color = (PERC_ALERTS_CORRECT_COLOR, 0.7),
            ),
            (
                char = :perc_true_outbreaks_detected,
                label = "Percent Outbreaks\nThat Are Detected",
                color = (PERC_OUTBREAKS_DETECTED_COLOR, 0.7)),
        ];
        columnfacetchar_label = "Detection Threshold",
        bins = -0.01:0.02:1.01,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_true_outbreak_alerts_perc_plot.png",
        ),
        compare_outbreak_true_outbreak_alerts_perc_plot;
        resolution = (2200, 1200),
    )

    compare_outbreak_cases_after_alerts_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
        (
            char = :cases_after_alerts,
            color = (PERC_OUTBREAKS_DETECTED_COLOR, 1.0),
        )
    ];
        xlabel = "Cases After Alerts",
        columnfacetchar_label = "Detection Threshold",
        binwidth = 50.0,
        meanlines = true,
        legend = false,
    )

    save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_cases_after_alerts_plot.png",
        ),
        compare_outbreak_cases_after_alerts_plot;
        resolution = (2200, 1200),
    )

    compare_outbreak_cases_perc_after_alerts_plot = compare_ensemble_OTchars_plots(
        ensemble_chars_vec,
        :detection_threshold,
        [
        (
            char = :cases_perc_after_alerts,
            color = (PERC_OUTBREAKS_DETECTED_COLOR, 1.0),
        )
    ];
        xlabel = "Percentage of Outbreak\nAfter Alerts",
        columnfacetchar_label = "Detection Threshold",
        binwidth = 0.01,
        meanlines = true,
        legend = false,
    )

    return save(
        plotsdir(
            "ensemble/testing-comparison/compare_outbreak_cases_perc_after_alerts_plot.png",
        ),
        compare_outbreak_cases_perc_after_alerts_plot;
        resolution = (2200, 1200),
    )
end
