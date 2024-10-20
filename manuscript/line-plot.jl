#%%
using DrWatson
@quickactivate "OutbreakDetection"

using OutbreakDetectionUtils
using OutbreakDetection
using Revise

include(srcdir("makie-plotting-setup.jl"))
include(srcdir("ensemble-parameters.jl"))

if false
    include("../src/ensemble-parameters.jl")
end

#%%
optimal_threshold_test_spec_vec = [
    IndividualTestSpecification(0.85, 0.85, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
    IndividualTestSpecification(1.0, 1.0, 14),
]

optimal_threshold_alertthreshold_vec = collect(1:1:15)

R_0_vec = [16.0]

ensemble_dynamics_spec_vec = create_combinations_vec(
    DynamicsParameters,
    (
        [ensemble_state_specification.init_states.N],
        [27],
        [0.2],
        [SIGMA],
        [GAMMA],
        R_0_vec,
        [0.8],
    ),
)

ensemble_spec_vec = create_combinations_vec(
    EnsembleSpecification,
    (
        [ensemble_model_type],
        [ensemble_state_specification],
        ensemble_dynamics_spec_vec,
        [ensemble_time_specification],
        [ensemble_nsims],
    ),
)

alert_method_vec = ["movingavg"]

#%%
ensemble_noise_specification =
    filter(ensemble_noise_specification_vec) do noise_spec
        noise_spec.noise_type == "poisson" ||
            noise_spec.correlation == "in-phase"
    end

ensemble_specification = ensemble_spec_vec[1]
alertmethod = alert_method_vec[1]

optimal_threshold_core_params = (
    alertthreshold_vec = optimal_threshold_alertthreshold_vec,
    ensemble_specification = ensemble_specification,
    outbreak_specification = ensemble_outbreak_specification,
    moving_avg_detection_lag = ensemble_moving_avg_detection_lag,
    percent_visit_clinic = ensemble_percent_visit_clinic,
    alertmethod = alertmethod,
)

basedirpath = joinpath(
    "R0_$(ensemble_specification.dynamics_parameters.R_0)"
    # noisespec_alertmethod_path,
)

baseplotdirpath = joinpath(
    plotsdir("ensemble/optimal-thresholds/lineplot"),
    basedirpath,
)

#%%
optimal_threshold_characteristics = collect_OptimalThresholdCharacteristics(
    ensemble_noise_specification,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params,
)

#%%
line_plot(
    optimal_threshold_characteristics;
    outcome = :accuracy,
    ylabel = "Outbreak Detection\nAccuracy",
    plotdirpath = baseplotdirpath,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0.5, 1.0),
    force = true,
    save_plot = false,
)

#%%
line_plot(
    optimal_threshold_characteristics;
    outcome = :unavoidable_cases,
    ylabel = "Unavoidable Cases",
    plotdirpath = baseplotdirpath,
    clinical_hline = false,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0, 3.5e4),
    force = true,
    save_plot = false,
)

#%%
line_plot(
    optimal_threshold_characteristics;
    outcome = :detectiondelays,
    ylabel = "Detection Delays\n(Days)",
    plotdirpath = baseplotdirpath,
    clinical_hline = false,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (-100, 100),
    force = true,
    save_plot = false,
)

#%%
line_plot(
    optimal_threshold_characteristics;
    outcome = :proportion_timeseries_in_outbreak,
    ylabel = "Proportion of Time Series\nIn Outbreak",
    plotdirpath = baseplotdirpath,
    clinical_hline = false,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0.0, 0.25),
    force = true,
    save_plot = false,
)

#%%
line_plot(
    optimal_threshold_characteristics;
    outcome = :proportion_timeseries_in_alert,
    ylabel = "Proportion of Time Series\nIn Alert",
    plotdirpath = baseplotdirpath,
    clinical_hline = false,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0.0, 0.35),
    force = true,
    save_plot = false,
)

#%%
line_plot(
    optimal_threshold_characteristics;
    outcome = :alert_outbreak_timeseries_prop_diff,
    ylabel = "Proportion of Time Series\nIn Alert - Outbreak",
    plotdirpath = baseplotdirpath,
    clinical_hline = false,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (-0.15, 0.20),
    force = true,
    save_plot = false,
)

#%%
accuracy_line_plot = line_plot(
    ensemble_noise_specification,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    outcome = :accuracy,
    ylabel = "Outbreak Detection\nAccuracy",
    plotdirpath = baseplotdirpath,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0.5, 1.0),
    force = true,
    save_plot = false,
)

#%%
unavoidable_line_plot = line_plot(
    ensemble_noise_specification,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    outcome = :unavoidable_cases,
    ylabel = "Unavoidable Cases",
    plotdirpath = baseplotdirpath,
    clinical_hline = false,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0, 3.5e4),
    force = true,
    save_plot = false,
)

#%%
delay_line_plot = line_plot(
    ensemble_noise_specification,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    outcome = :detectiondelays,
    ylabel = "Detection Delays\n(Days)",
    plotdirpath = baseplotdirpath,
    clinical_hline = false,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (-100, 100),
    force = true,
    save_plot = false,
)

#%%
outbreak_proportion_line_plot = line_plot(
    ensemble_noise_specification,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    outcome = :proportion_timeseries_in_outbreak,
    ylabel = "Proportion of Time Series\nIn Outbreak",
    plotdirpath = baseplotdirpath,
    clinical_hline = false,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (-1, 1),
    force = true,
    save_plot = false,
)

#%%
alert_proportion_line_plot = line_plot(
    ensemble_noise_specification,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    outcome = :proportion_timeseries_in_alert,
    ylabel = "Proportion of Time Series\nIn Alert",
    plotdirpath = baseplotdirpath,
    clinical_hline = false,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (-1, 1),
    force = true,
    save_plot = false,
)

#%%
alert_outbreak_proportion_line_plot = line_plot(
    ensemble_noise_specification,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    outcome = :alert_outbreak_timeseries_prop_diff,
    ylabel = "Proportion of Time Series\nIn Alert - Outbreak",
    plotdirpath = baseplotdirpath,
    clinical_hline = false,
    facet_fontsize = 18,
    labelsize = 20,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (-1, 1),
    force = true,
    save_plot = false,
)
