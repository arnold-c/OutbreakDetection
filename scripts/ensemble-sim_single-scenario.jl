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
ensemble_single_outbreak_detection_spec = OutbreakDetectionSpecification(
    5, 7, 0.6, 0.3
)

#%%
ensemble_single_seir_arr = get_ensemble_file(
    ensemble_specification
)["ensemble_seir_arr"]

ensemble_single_scenario_quantiles = get_ensemble_file(
    ensemble_specification, 95
)

ensemble_single_scenario_inc_file = get_ensemble_file(
    ensemble_specification, ensemble_outbreak_specification
)

ensemble_single_incarr = ensemble_single_scenario_inc_file["ensemble_inc_arr"]
ensemble_single_periodsum_vecs = ensemble_single_scenario_inc_file["ensemble_thresholds_vec"]

#%%
ensemble_single_scenario_quantiles_plot = create_sir_quantiles_plot(
    ensemble_single_scenario_quantiles["ensemble_seir_summary"];
    labels = seir_state_labels,
    colors = seircolors,
    annual = true,
    caption = ensemble_single_scenario_quantiles["caption"],
    timeparams = ensemble_time_specification,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble-sim_single-scenario_quantiles.png"
    ),
    ensemble_single_scenario_quantiles_plot,
)

#%%
ensemble_single_scenario_incidence_prevalence_plot = incidence_prevalence_plot(
    ensemble_single_incarr,
    ensemble_single_seir_arr,
    ensemble_single_periodsum_vecs,
    ensemble_time_specification;
    threshold = ensemble_outbreak_specification.outbreak_threshold,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_incidence_prevalence.png",
    ),
    ensemble_single_scenario_incidence_prevalence_plot,
)

#%%
for noise_specification in ensemble_noise_specification_vec
    scenario_specification = ScenarioSpecification(
        ensemble_specification,
        ensemble_outbreak_specification,
        noise_specification,
        ensemble_single_individual_test_spec,
    )

    ensemble_solution_dict = get_ensemble_file(
        scenario_specification.ensemble_specification
    )

    noisearr = create_noise_arr(
        noise_specification,
        ensemble_single_incarr;
        ensemble_specification = ensemble_specification,
        seed = 1234,
    )
    noisedir = getdirpath(noise_specification)

    plot_all_single_scenarios(
        noisearr,
        noisedir,
        OT_chars,
        ensemble_single_individual_test_spec,
        ensemble_single_outbreak_detection_spec,
        ensemble_time_specification,
    )
end
