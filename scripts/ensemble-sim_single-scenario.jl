#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using ColorSchemes

using OutbreakDetection

include(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
ensemble_single_individual_test_spec = IndividualTestSpecification(0.8, 0.8, 0)

ensemble_single_scenario_spec = ScenarioSpecification(
    ensemble_specification,
    ensemble_outbreak_specification,
    ensemble_noise_specification,
    OutbreakDetectionSpecification(5, 7, 0.6, 0.3),
    ensemble_single_individual_test_spec,
)

#%%
ensemble_single_scenario_sol = get_ensemble_file(
    ensemble_single_scenario_spec.ensemble_specification
)

ensemble_single_scenario_quantiles = get_ensemble_file(
    ensemble_single_scenario_spec.ensemble_specification, 95
)

ensemble_single_scenario_incarr = get_ensemble_file(
    ensemble_specification, ensemble_outbreak_specification
)

ensemble_single_scenario_detection = get_ensemble_file(
    ensemble_single_scenario_spec
)

#%%
ensemble_single_scenario_spec2 = ScenarioSpecification(
    ensemble_specification,
    ensemble_outbreak_specification,
    ensemble_noise_specification,
    OutbreakDetectionSpecification(10, 7, 0.6, 0.3),
    ensemble_single_individual_test_spec,
)

#%%
ensemble_single_scenario_detection2 = get_ensemble_file(
    ensemble_single_scenario_spec2
)

#%%
ensemble_single_scenario_noise_array = create_poisson_noise_arr(
    ensemble_single_scenario_incarr["ensemble_inc_arr"],
    ensemble_noise_specification,
)

#%%
ensemble_single_scenario_detection["testarr"] ==
ensemble_single_scenario_detection2["testarr"]

sum(ensemble_single_scenario_detection["testarr"][:, 7, 1])
sum(ensemble_single_scenario_detection2["testarr"][:, 7, 1])

ensemble_single_scenario_detection["OT_chars"].daily_sensitivity ==
ensemble_single_scenario_detection2["OT_chars"].daily_sensitivity

#%%
ensemble_single_scenario_quantiles_plot = create_sir_quantiles_plot(
    ensemble_single_scenario_quantiles["ensemble_seir_summary"];
    labels = seir_state_labels,
    colors = seircolors,
    annual = true,
    caption = ensemble_single_scenario_quantiles["caption"],
    timeparams = ensemble_single_scenario_spec.ensemble_specification.time_parameters,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_quantiles.png"
    ),
    ensemble_single_scenario_quantiles_plot,
)

#%%
ensemble_single_scenario_incidence_prevalence_plot = incidence_prevalence_plot(
    ensemble_single_scenario_incarr["ensemble_inc_arr"],
    ensemble_single_scenario_sol["ensemble_seir_arr"],
    ensemble_single_scenario_incarr["ensemble_thresholds_vec"],
    ensemble_single_scenario_spec.ensemble_specification.time_parameters;
    threshold = 5,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_incidence_prevalence.png",
    ),
    ensemble_single_scenario_incidence_prevalence_plot,
)

#%%
ensemble_single_scenario_noise_plot = visualize_ensemble_noise(
    ensemble_single_scenario_incarr["ensemble_inc_arr"],
    ensemble_noise_specification,
    ensemble_time_specification,
)

save(
    plotsdir("ensemble/single-scenario/ensemble-sim_single-scenario_noise.png"),
    ensemble_single_scenario_noise_plot,
)

#%%
ensemble_single_scenario_incidence_testing_plot = incidence_testing_plot(
    ensemble_single_scenario_incarr["ensemble_inc_arr"],
    ensemble_single_scenario_noise_array,
    ensemble_single_scenario_detection["testarr"],
    ensemble_single_scenario_spec.ensemble_specification.time_parameters,
    ensemble_single_scenario_spec.outbreak_detection_specification.alert_threshold;
    sim = 1,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_incidence-testing.png",
    ),
    ensemble_single_scenario_incidence_testing_plot,
)

#%%
ensemble_single_scenario_testing_timeseries_plot = testing_plot(
    ensemble_single_scenario_detection["testarr"],
    ensemble_single_scenario_spec.ensemble_specification.time_parameters,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_testing-timeseries.png",
    ),
    ensemble_single_scenario_testing_timeseries_plot,
)

#%%
ensemble_single_scenario_outbreak_dist_plot = ensemble_outbreak_distribution_plot(
    ensemble_single_scenario_detection["testarr"],
    ensemble_single_scenario_incarr["ensemble_inc_arr"],
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_outbreak-distribution.png",
    ),
    ensemble_single_scenario_outbreak_dist_plot,
)

#%%
ensemble_single_scenario_outbreak_alert_plot = ensemble_OTChars_plot(
    ensemble_single_scenario_detection["OT_chars"],
    ensemble_single_individual_test_spec,
    ensemble_single_scenario_spec.outbreak_detection_specification,
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
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_outbreak-alerts.png",
    ),
    ensemble_single_scenario_outbreak_alert_plot,
)

#%%
ensemble_single_scenario_outbreak_alert_perc_plot = ensemble_OTChars_plot(
    ensemble_single_scenario_detection["OT_chars"],
    ensemble_single_individual_test_spec,
    ensemble_single_scenario_spec.outbreak_detection_specification,
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
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_outbreak-alerts-perc.png",
    ),
    ensemble_single_scenario_outbreak_alert_perc_plot,
)

#%%
ensemble_single_scenario_outbreak_detect_diff_plot = ensemble_outbreak_detect_diff_plot(
    ensemble_single_scenario_detection["OT_chars"];
    binwidth = 1
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_outbreak-detect-diff.png",
    ),
    ensemble_single_scenario_outbreak_detect_diff_plot,
)

#%%
ensemble_single_scenario_sens_spec_dist_plot = ensemble_OTChars_plot(
    ensemble_single_scenario_detection["OT_chars"],
    ensemble_single_individual_test_spec,
    ensemble_single_scenario_spec.outbreak_detection_specification,
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
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_sens-spec-distribution.png",
    ),
    ensemble_single_scenario_sens_spec_dist_plot,
)

#%%
ensemble_single_scenario_ppv_npv_dist_plot = ensemble_OTChars_plot(
    ensemble_single_scenario_detection["OT_chars"],
    ensemble_single_individual_test_spec,
    ensemble_single_scenario_spec.outbreak_detection_specification,
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
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_ppv-npv-distribution.png",
    ),
    ensemble_single_scenario_ppv_npv_dist_plot,
)

#%%
ensemble_single_scenario_posodds_timeseries_plot = singlescenario_test_positivity_plot(
    ensemble_single_scenario_detection["test_positivity_structs"];
    agg = :thirty_day,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_posodds-timeseries.png",
    ),
    ensemble_single_scenario_posodds_timeseries_plot,
)

#%%
ensemble_single_scenario_posodds_dist_plot = test_positivity_distribution_plot(
    ensemble_single_scenario_detection["test_positivity_structs"];
    agg = :thirty_day,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_posodds-distribution.png",
    ),
    ensemble_single_scenario_posodds_dist_plot,
)

#%%
ensemble_single_scenario_posodds_outbreak_dist_plot = test_positivity_distribution_plot(
    ensemble_single_scenario_detection["test_positivity_structs"];
    agg = :seven_day,
    color = :outbreak => "Outbreak Status",
    layout = :outbreak,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_posodds-outbreak-distribution.png",
    ),
    ensemble_single_scenario_posodds_outbreak_dist_plot,
)

#%%
ensemble_single_scenario_detection_delay_dist_plot = ensemble_OTChars_plot(
    ensemble_single_scenario_detection["OT_chars"],
    ensemble_single_individual_test_spec,
    ensemble_single_scenario_spec.outbreak_detection_specification,
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
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_detection-delay-distribution.png",
    ),
    ensemble_single_scenario_detection_delay_dist_plot,
)
