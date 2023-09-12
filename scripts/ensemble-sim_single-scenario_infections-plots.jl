#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ColorSchemes

include("ensemble-sim_single-scenario.jl")

include(srcdir("makie-plotting-setup.jl"))

#%%
create_sir_quantiles_plot(
    ensemble_single_scenario_quantiles["ensemble_seir_summary"];
    labels = seir_state_labels,
    colors = seircolors,
    annual = true,
    caption = caption,
    timeparams = ensemble_single_scenario_spec.ensemble_specification.time_parameters,
)

#%%
outbreakcols = [ColorSchemes.magma[i] for i in (200, 20)]

detect_outbreak_plot(
    ensemble_single_scenario_detection["incarr"],
    ensemble_single_scenario_sol["ensemble_seir_arr"],
    ensemble_single_scenario_spec.ensemble_specification.time_parameters;
    colormap = outbreakcols,
    # xlims = (90, 100),
    # ylims_inc = (0, 150),
    # ylims_periodsum = (0, 1000),
)

#%%
visualize_ensemble_noise(
    ensemble_single_scenario_spec.noise_specification.noise_array,
    ensemble_single_scenario_spec.noise_specification.time_parameters,
)

#%%
incidence_testing_plot(
    ensemble_single_scenario_detection["incarr"],
    ensemble_single_scenario_detection["testarr"],
    ensemble_single_scenario_spec.ensemble_specification.time_parameters,
    ensemble_single_scenario_spec.outbreak_detection_specification.detection_threshold;
    sim = 1,
)

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
