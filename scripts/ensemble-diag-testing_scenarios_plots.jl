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
detectthreshold_vec = [collect(5:15:20)..., 35, 50, 100]

#%%
ensemble_chars_vec = Vector(
    undef, length(sensitivity_vec) * length(detectthreshold_vec)
)

#%%
ensemble_static_noise_arr = create_static_NoiseSpecification(
    [10.0],
    SimTimeParameters(;
        tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
    ),
    0.0,
    0.1,
    1_000,
)

#%%
prog = Progress(length(sensitivity_vec) * length(detectthreshold_vec))
@floop for (i, ((sens, spec), detectthrehold)) in enumerate(
    Iterators.product(
        zip(sensitivity_vec, specificity_vec),
        detectthreshold_vec
    ),
)
    ind_test_spec = IndividualTestSpecification(sens, spec)
    outbreak_detect_spec = OutbreakDetectionSpecification(
        detectthrehold, 7, 0.3, 0.3, 3
    )

    ensemble_spec =
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
                    DynamicsParameters(500_000, 10, 0.2),
                    time_params,
                    1_000,
                ),
                OutbreakSpecification(5, 30, 500),
                ensemble_static_noise_arr,
                outbreak_detect_spec,
                ind_test_spec,
            )
        end

    ensemble_chars_file = OutbreakDetection.get_scenario_file("scenario", ensemble_spec)

    ensemble_chars_vec[i] = (
        OT_chars = ensemble_chars_file["OT_chars"],
        outbreak_detect_spec = outbreak_detect_spec,
        ind_test_spec = ind_test_spec,
    )

    next!(prog)
end

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

save(plotsdir("ensemble/testing-comparison/compare_outbreak_sens_spec_plot.png"), compare_outbreak_sens_spec_plot)

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
    char2_color = :purple
)

save(plotsdir("ensemble/testing-comparison/compare_outbreak_ppv_npv_plot.png"), compare_outbreak_ppv_npv_plot)
