#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using FLoops
using NaNMath: NaNMath

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
    :detection_threshold;
    columnfacetchar_label = "Detection Threshold",
    bins = 0.0:0.01:1.01,
    plottingchars = [
        (
            char = :daily_sensitivity,
            label = "Sensitivity",
            color = (:blue, 0.5),
        ),
        (
            char = :daily_specificity,
            label = "Specificity",
            color = (:red, 0.5),
        ),
    ],
)

save(
    plotsdir("ensemble/testing-comparison/compare_outbreak_sens_spec_plot.png"),
    compare_outbreak_sens_spec_plot;
    resolution = (2200, 1200),
);

#%%
compare_outbreak_ppv_npv_plot = compare_ensemble_OTchars_plots(
    ensemble_chars_vec,
    :detection_threshold;
    columnfacetchar_label = "Detection Threshold",
    bins = 0.0:0.01:1.01,
    plottingchars = [
        (
            char = :daily_ppv,
            label = "PPV",
            color = (:green, 0.5),
        ),
        (
            char = :daily_npv,
            label = "NPV",
            color = (:purple, 0.5),
        ),
    ],
)

save(
    plotsdir("ensemble/testing-comparison/compare_outbreak_ppv_npv_plot.png"),
    compare_outbreak_ppv_npv_plot;
    resolution = (2200, 1200),
);

#%%
compare_outbreak_detection_delays_plot = compare_ensemble_OTchars_plots(
    ensemble_chars_vec,
    :detection_threshold;
    columnfacetchar_label = "Detection Threshold",
    binwidth = 5.0,
    plottingchars = [
    (
        char = :detectiondelays,
        label = "Detection Delay",
        color = (:red, 0.8),
    )
    ],
)

save(
    plotsdir(
        "ensemble/testing-comparison/compare_outbreak_detection_delays_plot.png"
    ),
    compare_outbreak_detection_delays_plot;
    resolution = (2200, 1200),
);

#%%
compare_outbreak_alert_numbers_plot = compare_ensemble_OTchars_plots(
    ensemble_chars_vec,
    :detection_threshold;
    columnfacetchar_label = "Detection Threshold",
    binwidth = 1.0,
    plottingchars = [
    (
        char = :n_alerts_per_outbreak,
        label = "Alerts per Outbreak",
        color = (:orange, 1.0),
    )
    ],
)

save(
    plotsdir(
        "ensemble/testing-comparison/compare_outbreak_alert_numbers_plot.png"
    ),
    compare_outbreak_alert_numbers_plot;
    resolution = (2200, 1200),
);

#%%
compare_outbreak_false_alerts_plot = compare_ensemble_OTchars_plots(
    ensemble_chars_vec,
    :detection_threshold;
    columnfacetchar_label = "Detection Threshold",
    binwidth = 1.0,
    plottingchars = [
    (
        char = :n_false_alerts,
        label = "False Alerts",
        color = (:navy, 1.0),
    )
    ],
)

save(
    plotsdir(
        "ensemble/testing-comparison/compare_outbreak_false_alerts_plot.png"
    ),
    compare_outbreak_false_alerts_plot;
    resolution = (2200, 1200),
)

#%%
compare_outbreak_missed_plot = compare_ensemble_OTchars_plots(
    ensemble_chars_vec,
    :detection_threshold;
    columnfacetchar_label = "Detection Threshold",
    binwidth = 1.0,
    plottingchars = [
    (
        char = :n_missed_outbreaks,
        label = "Missed Outbreaks",
        color = (:grey20, 1.0),
    )
    ],
)

save(
    plotsdir(
        "ensemble/testing-comparison/compare_outbreak_missed_plot.png"
    ),
    compare_outbreak_missed_plot;
    resolution = (2200, 1200),
)

#%%
compare_outbreak_true_outbreak_perc_plot = compare_ensemble_OTchars_plots(
    ensemble_chars_vec,
    :detection_threshold;
    columnfacetchar_label = "Detection Threshold",
    binwidth = 0.02,
    plottingchars = [
        (
            char = :perc_true_outbreaks_detected,
            label = "Percent True Outbreaks Detected",
            color = (:navy, 0.5),
        ),
        (
            char = :perc_true_outbreaks_missed,
            label = "Percent True Outbreaks Missed",
            color = (:orange, 0.5)),
    ],
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

#%%
compare_outbreak_alerts_perc_plot = compare_ensemble_OTchars_plots(
    ensemble_chars_vec,
    :detection_threshold;
    columnfacetchar_label = "Detection Threshold",
    bins = -0.01:0.02:1.01,
    plottingchars = [
        (
            char = :perc_alerts_correct,
            label = "Percent Alerts\nThat Are Correct",
            color = (:green, 0.5),
        ),
        (
            char = :perc_alerts_false,
            label = "Percent Alerts\nThat Are False",
            color = (:grey20, 0.5)),
    ],
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

#%%
compare_outbreak_true_outbreak_alerts_perc_plot = compare_ensemble_OTchars_plots(
    ensemble_chars_vec,
    :detection_threshold;
    columnfacetchar_label = "Detection Threshold",
    bins = -0.01:0.02:1.01,
    plottingchars = [
        (
            char = :perc_alerts_correct,
            label = "Percent Alerts\nThat Are Correct",
            color = (:green, 0.5),
        ),
        (
            char = :perc_true_outbreaks_detected,
            label = "Percent True Outbreaks\nThat Are Detected",
            color = (:navy, 0.5)),
    ],
)

save(
    plotsdir(
        "ensemble/testing-comparison/compare_outbreak_true_outbreak_alerts_perc_plot.png",
    ),
    compare_outbreak_true_outbreak_alerts_perc_plot;
    resolution = (2200, 1200),
)
