#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using FLoops

using OutbreakDetection

include(srcdir("makie-plotting-setup.jl"))

#%%
sensitivity_vec = collect(0.8:0.2:1.0)
specificity_vec = collect(0.8:0.2:1.0)
detectthreshold_vec = [2, 4, collect(5:5:20)...]

#%%
ensemble_scenario_spec_vec = Vector{ScenarioSpecification}(
    undef,
    length(sensitivity_vec) * length(detectthreshold_vec) +
    length(detectthreshold_vec),
)
ensemble_chars_vec = Vector(
    undef, length(ensemble_scenario_spec_vec)
)

#%%
ensemble_specification = EnsembleSpecification(
    ("seasonal-infectivity-import", "tau-leaping"),
    StateParameters(
        500_000,
        Dict(
            :s_prop => 0.1,
            :e_prop => 0.0,
            :i_prop => 0.0,
            :r_prop => 0.9,
        ),
    ),
    DynamicsParameters(500_000, 10, 0.2),
    SimTimeParameters(;
        tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
    ),
    100,
)
noise_specification = NoiseSpecification("poisson", 1.0)
outbreak_specification = OutbreakSpecification(5, 30, 500)

percent_visit_clinic = 0.6
outbreak_detect_spec_vec = map(
    threshold -> OutbreakDetectionSpecification(
        threshold, 7, percent_visit_clinic, 0.8, 0
    ),
    detectthreshold_vec,
)
clinical_case_outbreak_detect_spec_vec = map(
    threshold -> OutbreakDetectionSpecification(
        threshold, 7, percent_visit_clinic, 1.0, 0
    ),
    detectthreshold_vec,
)

#%%
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

#%%
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

#%%
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

#%%
sort!(
    ensemble_chars_vec;
    by = x -> (
        x.outbreak_detect_spec.detection_threshold,
        x.ind_test_spec.specificity,
    ),
);

#%%
compare_outbreak_sens_spec_plot = compare_ensemble_OTchars_plots(
    ensemble_chars_vec,
    :sensitivity,
    :specificity,
    :detection_threshold;
    char1_label = "Sensitivity",
    char2_label = "Specificity",
    char3_label = "Detection Threshold",
)

save(
    plotsdir("ensemble/testing-comparison/compare_outbreak_sens_spec_plot.png"),
    compare_outbreak_sens_spec_plot;
    resolution = (2200, 1200),
);

#%%
compare_outbreak_ppv_npv_plot = compare_ensemble_OTchars_plots(
    ensemble_chars_vec,
    :ppv,
    :npv,
    :detection_threshold;
    char1_label = "PPV",
    char2_label = "NPV",
    char3_label = "Detection Threshold",
    char1_color = :green,
    char2_color = :purple,
)

save(
    plotsdir("ensemble/testing-comparison/compare_outbreak_ppv_npv_plot.png"),
    compare_outbreak_ppv_npv_plot;
    resolution = (2200, 1200),
);

#%%
compare_outbreak_detection_delays_plot = compare_ensemble_OTchars_plots(
    ensemble_chars_vec,
    :detectiondelays,
    :detection_threshold;
    char1_label = "Detection Delay",
    char2_label = "Detection Threshold",
    char1_color = :blue,
    bins = -11.0:1.0:150.0,
)

save(
    plotsdir(
        "ensemble/testing-comparison/compare_outbreak_detection_delays_plot.png"
    ),
    compare_outbreak_detection_delays_plot;
    resolution = (2200, 1200),
);
