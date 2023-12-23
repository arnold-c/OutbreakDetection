function plot_all_single_scenarios(
    noisearr,
    noisedir,
    OT_chars,
    incarr,
    testarr,
    test_specification,
    outbreak_detection_specification,
    time_specification,
)
    ensemble_noise_plotpath = joinpath(
        plotsdir(),
        "ensemble",
        "single-scenario",
        noisedir
    )
    mkpath(ensemble_noise_plotpath)

    ensemble_noise_fig = visualize_ensemble_noise(
        noisearr,
        time_specification,
        noisedir
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_noise.png"
        ),
        ensemble_noise_fig; resolution = (2200, 1600),
    )

    ensemble_single_scenario_outbreak_alert_plot = ensemble_OTChars_plot(
        OT_chars,
        test_specification,
        outbreak_detection_specification,
        (
            (
                char = :noutbreaks,
                label = "Number of Outbreaks",
                color = (N_OUTBREAKS_COLOR, 0.5),
                hjust = 11,
                vjust = 50,
            ),
            (
                char = :nalerts,
                label = "Number of Alerts",
                color = (N_ALERTS_COLOR, 0.5),
                hjust = 11,
                vjust = 50,
            ),
        ),
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_outbreak-alerts.png",
        ),
        ensemble_single_scenario_outbreak_alert_plot,
    )

    ensemble_single_scenario_outbreak_alert_perc_plot = ensemble_OTChars_plot(
        OT_chars,
        test_specification,
        outbreak_detection_specification,
        (
            (
                char = :perc_true_outbreaks_detected,
                label = "Percent Outbreaks Detected",
                color = (PERC_OUTBREAKS_DETECTED_COLOR, 1.0),
                hjust = -0.15,
                vjust = 60,
            ),
            (
                char = :perc_alerts_correct,
                label = "Percent Alerts\nThat Are Correct",
                color = (PERC_ALERTS_CORRECT_COLOR, 0.7),
                hjust = -0.1,
                vjust = 60,
            ),
        );
        bins = 0.0:0.01:1.01,
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_outbreak-alerts-perc.png",
        ),
        ensemble_single_scenario_outbreak_alert_perc_plot,
    )

    #%%
    ensemble_single_scenario_outbreak_detect_diff_plot = ensemble_outbreak_detect_diff_plot(
        OT_chars;
        binwidth = 1
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_outbreak-detect-diff.png",
        ),
        ensemble_single_scenario_outbreak_detect_diff_plot,
    )

    ensemble_single_scenario_sens_spec_dist_plot = ensemble_OTChars_plot(
        OT_chars,
        test_specification,
        outbreak_detection_specification,
        (
            (
                char = :daily_sensitivity,
                label = "Sensitivity",
                color = (DAILY_SENSITIVITY_COLOR, 0.5),
                hjust = -0.085,
                vjust = 80,
            ),
            (
                char = :daily_specificity,
                label = "Specificity",
                color = (DAILY_SPECIFICITY_COLOR, 0.5),
                hjust = -0.08,
                vjust = 75,
            ),
        );
        bins = -0.005:0.01:1.005,
        legendlabel = "Outbreak Characteristic",
        normalization = :pdf,
        meanlines = true,
        meanlabels = true,
        meanannotations = true,
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_sens-spec-distribution.png",
        ),
        ensemble_single_scenario_sens_spec_dist_plot,
    )

    ensemble_single_scenario_ppv_npv_dist_plot = ensemble_OTChars_plot(
        OT_chars,
        test_specification,
        outbreak_detection_specification,
        (
            (
                char = :daily_ppv,
                label = "PPV",
                color = (DAILY_PPV_COLOR, 0.5),
                hjust = 0.01,
                vjust = 80,
            ),
            (
                char = :daily_npv,
                label = "NPV",
                color = (DAILY_NPV_COLOR, 0.5),
                hjust = -0.05,
                vjust = 75,
            ),
        );
        bins = -0.005:0.01:1.005,
        normalization = :pdf,
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_ppv-npv-distribution.png",
        ),
        ensemble_single_scenario_ppv_npv_dist_plot,
    )

    ensemble_single_scenario_detection_delay_dist_plot = ensemble_OTChars_plot(
        OT_chars,
        test_specification,
        outbreak_detection_specification,
        (
        (
            char = :detectiondelays,
            label = "Detection Delay",
            color = (DETECTION_DELAY_COLOR, 0.8),
            hjust = -50,
            vjust = 700,
        ),
    );
        binwidth = 1,
        xlabel = "Detection Delay (days)",
        legend = false,
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_detection-delay-distribution.png",
        ),
        ensemble_single_scenario_detection_delay_dist_plot,
    )

    incidence_testing_plottitle = "Sens: $(test_specification.sensitivity), Spec: $(test_specification.specificity), Lag: $(test_specification.test_result_lag),\nThreshold: $(outbreak_detection_specification.alert_threshold), Perc Clinic Tested: $(outbreak_detection_specification.percent_clinic_tested)\nNoise: $(noisedir)"

    ensemble_single_scenario_incidence_testing_plot = incidence_testing_plot(
        incarr,
        noisearr,
        testarr,
        outbreak_detection_specification,
        time_specification;
        sim = 1,
        plottitle = incidence_testing_plottitle,
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_incidence-testing.png",
        ),
        ensemble_single_scenario_incidence_testing_plot,
    )

    ensemble_single_scenario_testing_timeseries_plot = testing_plot(
        testarr,
        time_specification;
        plottitle = incidence_testing_plottitle
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_testing-timeseries.png",
        ),
        ensemble_single_scenario_testing_timeseries_plot,
    )

    # ensemble_single_scenario_outbreak_dist_plot = ensemble_outbreak_distribution_plot(
    #     ensemble_single_scenario_detection["testarr"],
    #     ensemble_single_scenario_incarr["ensemble_inc_arr"],
    # )
    #
    # save(
    #     plotsdir(
    #         "ensemble/single-scenario/ensemble-sim_single-scenario_outbreak-distribution.png",
    #     ),
    #     ensemble_single_scenario_outbreak_dist_plot,
    # )

    return nothing
end
