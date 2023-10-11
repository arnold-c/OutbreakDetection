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
    DynamicsParameters(500_000, 10, 0.2; vaccination_coverage = 0.8),
    ensemble_single_time_spec,
    10,
)

ensemble_single_outbreak_spec = OutbreakSpecification(5, 30, 500)
ensemble_single_individual_test_spec = IndividualTestSpecification(0.8, 0.8)
ensemble_single_noise_spec = NoiseSpecification("poisson", 1.0)
ensemble_single_scenario_spec = ScenarioSpecification(
    ensemble_single_ensemble_spec,
    ensemble_single_outbreak_spec,
    ensemble_single_noise_spec,
    OutbreakDetectionSpecification(4, 7, 0.6, 0.8, 3),
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
    OutbreakDetectionSpecification(10, 7, 0.6, 0.8, 3),
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

ensemble_single_scenario_detection["OT_chars"].sensitivity ==
ensemble_single_scenario_detection2["OT_chars"].sensitivity

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
ensemble_single_scenario_detect_outbreak_plot = detect_outbreak_plot(
    ensemble_single_scenario_incarr["ensemble_inc_arr"],
    ensemble_single_scenario_sol["ensemble_seir_arr"],
    ensemble_single_scenario_spec.ensemble_specification.time_parameters;
    colormap = outbreakcols,
    # xlims = (90, 100),
    # ylims_inc = (0, 150),
    # ylims_periodsum = (0, 1000),
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_detect-outbreak.png",
    ),
    ensemble_single_scenario_detect_outbreak_plot,
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
    :noutbreaks,
    :ndetectoutbreaks,
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
    :sensitivity,
    :specificity;
    bins = 0.0:0.01:1.01,
    char1_label = "Sensitivity",
    char2_label = "Specificity",
    char1_color = :red,
    char2_color = :blue,
    xlabel = "Characteristic Value",
    legendlabel = "Outbreak Characteristic",
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
    :ppv,
    :npv;
    bins = 0.0:0.01:1.01,
    char1_label = "PPV",
    char2_label = "NPV",
    char1_color = :green,
    char2_color = :purple,
    xlabel = "Characteristic Value",
    legendlabel = "Outbreak Characteristic",
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_ppv-npv-distribution.png",
    ),
    ensemble_single_scenario_ppv_npv_dist_plot,
)

#%%
ensemble_single_scenario_posodds_timeseries_plot = singlescenario_test_positivity_plot(
    [ensemble_single_scenario_detection["test_positivity_structs"][1]];
    agg = :thirty_day
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_posodds-timeseries.png",
    ),
    ensemble_single_scenario_posodds_timeseries_plot,
)

#%%
ensemble_single_scenario_posodds_dist_plot = @chain ensemble_single_scenario_detection["posoddsarr"][:, 1, :] begin
    replace(_, NaN => 0.0)
    test_positivity_distribution_plot(_)
end

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_posodds-distribution.png",
    ),
    ensemble_single_scenario_posodds_dist_plot,
)
