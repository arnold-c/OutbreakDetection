#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using ColorSchemes

using OutbreakDetection

include(srcdir("makie-plotting-setup.jl"))

#%%
ensemble_single_scenario_spec =
    let time_params = SimTimeParameters(;
            tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
        )
        ScenarioSpecification(
            EnsembleSpecification(
                ("seasonal-infectivity-import", "tau-leaping"),
                StateParameters(
                    500_000,
                    Dict(
                        :s_prop => 0.1,
                        :e_prop => 0.01,
                        :i_prop => 0.01,
                        :r_prop => 0.88,
                    ),
                ),
                DynamicsParameters(500_000, 5, 0.2),
                time_params,
                1_000,
            ),
            OutbreakSpecification(5, 30, 500),
            create_static_NoiseSpecification(
                [10.0],
                time_params,
                0.0,
                0.1,
                1_000
            ),
            OutbreakDetectionSpecification(10, 7, 0.3, 0.3, 3),
            IndividualTestSpecification(0.8, 0.8),
        )
    end

#%%
ensemble_single_scenario_sol = get_ensemble_file(
    "solution", ensemble_single_scenario_spec.ensemble_specification
)

ensemble_single_scenario_quantiles = get_ensemble_file(
    "95", ensemble_single_scenario_spec.ensemble_specification
)

ensemble_single_scenario_detection = get_scenario_file(
    "scenario", ensemble_single_scenario_spec
)

#%%
ensemble_single_scenario_quantiles_plot = create_sir_quantiles_plot(
    ensemble_single_scenario_quantiles["ensemble_seir_summary"];
    labels = seir_state_labels,
    colors = seircolors,
    annual = true,
    caption = ensemble_single_scenario_quantiles["caption"],
    timeparams = ensemble_single_scenario_spec.ensemble_specification.time_parameters,
)

save(plotsdir("ensemble-sim_single-scenario_quantiles.png"), ensemble_single_scenario_quantiles_plot)

#%%
outbreakcols = [ColorSchemes.magma[i] for i in (200, 20)]

ensemble_single_scenario_detect_outbreak_plot = detect_outbreak_plot(
    ensemble_single_scenario_detection["incarr"],
    ensemble_single_scenario_sol["ensemble_seir_arr"],
    ensemble_single_scenario_spec.ensemble_specification.time_parameters;
    colormap = outbreakcols,
    # xlims = (90, 100),
    # ylims_inc = (0, 150),
    # ylims_periodsum = (0, 1000),
)

save(plotsdir("ensemble-sim_single-scenario_detect-outbreak.png"), ensemble_single_scenario_detect_outbreak_plot)

#%%
ensemble_single_scenario_noise_plot = visualize_ensemble_noise(
    ensemble_single_scenario_spec.noise_specification.noise_array,
    ensemble_single_scenario_spec.noise_specification.time_parameters,
)

save(plotsdir("ensemble-sim_single-scenario_noise.png"), ensemble_single_scenario_noise_plot)

#%%
ensemble_single_scenario_incidence_testing_plot = incidence_testing_plot(
    ensemble_single_scenario_detection["incarr"],
    ensemble_single_scenario_detection["testarr"],
    ensemble_single_scenario_spec.ensemble_specification.time_parameters,
    ensemble_single_scenario_spec.outbreak_detection_specification.detection_threshold;
    sim = 1,
)

save(plotsdir("ensemble-sim_single-scenario_incidence-testing.png"), ensemble_single_scenario_incidence_testing_plot)

#%%
testing_plot(
    ensemble_single_scenario_detection["testarr"],
    ensemble_single_scenario_spec.ensemble_specification.time_parameters,
)

#%%
ensemble_outbreak_distribution_plot(
    ensemble_single_scenario_detection["testarr"],
    ensemble_single_scenario_detection["incarr"],
)

#%%
ensemble_OTChars_plot(
    ensemble_single_scenario_detection["OT_chars"],
    :sensitivity,
    :specificity;
    bins = 0.0:0.01:1.01,
    char1_label = "Sensitivity",
    char2_label = "Specificity",
    xlabel = "Proportion",
    legendlabel = "Characteristic",
)

#%%
ensemble_OTChars_plot(
    ensemble_single_scenario_detection["OT_chars"],
    :noutbreaks,
    :ndetectoutbreaks,
)

#%%
ensemble_OTChars_plot(
    ensemble_single_scenario_detection["OT_chars"],
    :ppv,
    :npv;
    bins = 0.0:0.01:1.01,
    char1_label = "PPV",
    char2_label = "NPV",
    char1_color = :green,
    char2_color = :purple,
    xlabel = "Proportion",
    legendlabel = "Characteristic",
)
