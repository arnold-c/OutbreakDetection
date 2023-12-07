#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using FLoops
using NaNMath: NaNMath
using DataFrames
using DataFramesMeta
using Statistics

using OutbreakDetection

includet(srcdir("makie-plotting-setup.jl"))

#%%
test_spec_vec = [
    IndividualTestSpecification(0.5, 0.5, 0),
    IndividualTestSpecification(0.7, 0.7, 0),
    IndividualTestSpecification(0.8, 0.8, 0),
    IndividualTestSpecification(0.85, 0.85, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    CLINICAL_TEST_SPECS...,
    IndividualTestSpecification(1.0, 1.0, 0),
    IndividualTestSpecification(1.0, 1.0, 3),
    IndividualTestSpecification(1.0, 1.0, 7),
    IndividualTestSpecification(1.0, 1.0, 14),
]

alertthreshold_vec = collect(4:1:30)

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
    DynamicsParameters(500_000, 10, 0.2; vaccination_coverage = 0.0),
    SimTimeParameters(;
        tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
    ),
    100,
)
noise_specification = NoiseSpecification("poisson", 1.0)
outbreak_specification = OutbreakSpecification(5, 30, 500)

moving_avg_detection_lag = 7
percent_visit_clinic = 0.6
percent_clinic_tested_vec = collect(0.2:0.2:1.0)

threshold_comparison_params = (
    alertthreshold_vec = alertthreshold_vec,
    ensemble_specification = ensemble_specification,
    noise_specification = noise_specification,
    outbreak_specification = outbreak_specification,
    moving_avg_detection_lag = moving_avg_detection_lag,
    percent_visit_clinic = percent_visit_clinic,
)

#%%
optimal_thresholds_vec = calculate_OptimalThresholdCharacteristics(
    percent_clinic_tested_vec,
    test_spec_vec,
    threshold_comparison_params
)

#%%
compare_optimal_thresholds_chars_plot(
    optimal_thresholds_vec,
    [
        (
            char = :accuracy,
            label = "Accuracy",
            color = (ACCURACY_COLOR, 0.7),
            binwidth = 0.01,
        ),
        (
            char = :detectiondelays,
            label = "Detection Delay (Days)",
            color = (DETECTION_DELAY_COLOR, 1.0),
            binwidth = 10,
        ),
        (
            char = :cases_before_alerts,
            label = "Number of Cases Before Alert",
            color = (PERC_OUTBREAKS_MISSED_COLOR, 1.0),
            binwidth = 50,
        ),
        (
            char = :cases_after_alerts,
            label = "Number of Cases After Alert",
            color = (PERC_OUTBREAKS_DETECTED_COLOR, 1.0),
            binwidth = 500,
        ),
        (
            char = :detected_outbreak_size,
            label = "Size of Outbreaks Detected",
            color = (PERC_OUTBREAKS_DETECTED_COLOR, 0.7),
            binwidth = 500,
        ),
        (
            char = :missed_outbreak_size,
            label = "Size of Outbreaks Missed",
            color = (PERC_OUTBREAKS_MISSED_COLOR, 0.7),
            binwidth = 100,
        ),
    ],
)

#%%
create_and_save_xlsx_optimal_threshold_summaries(optimal_thresholds_vec)

create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec, :detectiondelays
)

create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec, :cases_before_alerts
)

create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec, :cases_after_alerts
)

create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec, :detected_outbreak_size
)

create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec, :missed_outbreak_size
)
