#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using ColorSchemes

using OutbreakDetection

include(srcdir("makie-plotting-setup.jl"))

#%%
ensemble_single_time_spec = SimTimeParameters(;
    tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
)

ensemble_single_ensemble_spec = EnsembleSpecification(
    ("seasonal-infectivity-import", "tau-leaping"),
    StateParameters(
        500_000,
        Dict(:s_prop => 0.1, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.9),
    ),
    DynamicsParameters(500_000, 10, 0.2; vaccination_coverage = 0.0),
    ensemble_single_time_spec,
    100,
)

ensemble_single_outbreak_spec = OutbreakSpecification(5, 30, 500)
ensemble_single_individual_test_spec = IndividualTestSpecification(0.8, 0.8)
ensemble_single_noise_spec = NoiseSpecification("poisson", 1.0)
ensemble_single_scenario_spec = ScenarioSpecification(
    ensemble_single_ensemble_spec,
    ensemble_single_outbreak_spec,
    ensemble_single_noise_spec,
    OutbreakDetectionSpecification(5, 7, 0.6, 0.8, 0),
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
    ensemble_single_ensemble_spec, ensemble_single_outbreak_spec
)

ensemble_single_scenario_detection = get_ensemble_file(
    ensemble_single_scenario_spec
)

#%%
ensemble_single_scenario_spec2 = ScenarioSpecification(
    ensemble_single_ensemble_spec,
    ensemble_single_outbreak_spec,
    ensemble_single_noise_spec,
    OutbreakDetectionSpecification(10, 7, 0.6, 0.8, 0),
    ensemble_single_individual_test_spec,
)

#%%
ensemble_single_scenario_detection2 = get_ensemble_file(
    ensemble_single_scenario_spec2
)

#%%
ensemble_single_scenario_noise_array = create_poisson_noise_arr(
    ensemble_single_scenario_incarr["ensemble_inc_arr"],
    ensemble_single_noise_spec,
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
    ensemble_single_noise_spec,
    ensemble_single_time_spec,
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
    ensemble_single_scenario_spec.outbreak_detection_specification.detection_threshold;
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
ensemble_single_scenario_outbreak_detect_plot = ensemble_OTChars_plot(
    ensemble_single_scenario_detection["OT_chars"],
    ensemble_single_individual_test_spec,
    ensemble_single_scenario_spec.outbreak_detection_specification,
    (
        (
            char = :noutbreaks,
            label = "Number of Outbreaks",
            color = (N_OUTBREAKS_COLOR, 0.5),
            hjust = 11,
            vjust = 0.055,
        ),
        (
            char = :ndetectoutbreaks,
            label = "Number of Alerts",
            color = (N_ALERTS_COLOR, 0.5),
            hjust = 2.5,
            vjust = 0.055,
        ),
    ),
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_outbreak-detection.png",
    ),
    ensemble_single_scenario_outbreak_detect_plot,
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
            hjust = -0.085,
            vjust = 80,
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
            hjust = -0.07,
            vjust = 80,
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
        hjust = -18,
        vjust = 0.18,
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
