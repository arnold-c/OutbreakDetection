#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using ColorSchemes

using OutbreakDetection

includet(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
ensemble_single_individual_test_spec = IndividualTestSpecification(
    0.85, 0.85, 0
)
ensemble_single_outbreak_detection_spec = OutbreakDetectionSpecification(
    5, 7, 0.6, 0.3, "dailythreshold_movingavg"
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

mkpath(plotsdir("ensemble/single-scenario"))

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
create_testing_related_plots(
    ensemble_specification,
    ensemble_outbreak_specification,
    ensemble_noise_specification_vec[2],
    ensemble_single_outbreak_detection_spec,
    ensemble_single_individual_test_spec,
    ensemble_time_specification,
    ensemble_single_incarr;
    seed = 1234,
    sim = 1,
    xlims = (0, 10),
)

#%%
for noise_specification in ensemble_noise_specification_vec
    create_testing_related_plots(
        ensemble_specification,
        ensemble_outbreak_specification,
        noise_specification,
        ensemble_single_outbreak_detection_spec,
        ensemble_single_individual_test_spec,
        ensemble_time_specification,
        ensemble_single_incarr;
        seed = 1234,
        sim = 1,
    )

    GC.gc(true)
    @info "Finished plotting the single scenario for $(noisedir)"
    println("=================================================================")
end
