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
    IndividualTestSpecification(0.8, 0.8, 0),
    CLINICAL_CASE_TEST_SPEC,
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
optimal_thresholds_df = DataFrame(;
    percent_clinic_tested = optimal_thresholds_vec.percent_clinic_tested,
    sensitivity = getfield.(
        optimal_thresholds_vec.individual_test_specification, :sensitivity
    ),
    specificity = getfield.(
        optimal_thresholds_vec.individual_test_specification, :specificity
    ),
    alert_threshold = optimal_thresholds_vec.alert_threshold,
    accuracy = optimal_thresholds_vec.accuracy,
)

#%%
@chain optimal_thresholds_df begin
    @orderby :specificity
    unstack(
        _,
        [:sensitivity, :specificity, :test_lag],
        :percent_clinic_tested,
        :alert_threshold,
    )
    select(_, Cols(x -> startswith(x, "s"), x -> startswith(x, "0"), "1.0"))
end

#%%
@chain optimal_thresholds_df begin
    @orderby :specificity
    unstack(
        _,
        [:sensitivity, :specificity, :test_lag],
        :percent_clinic_tested,
        :accuracy,
    )
    select(_, Cols(x -> startswith(x, "s"), x -> startswith(x, "0"), "1.0"))
end

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
            label = "Number of Cases After Alert",
            color = (PERC_OUTBREAKS_MISSED_COLOR, 1.0),
            binwidth = 500,
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
detectiondelays_df = create_optimal_threshold_summary_df(:detectiondelays)
cases_before_alerts_df = create_optimal_threshold_summary_df(
    :cases_before_alerts
)
cases_after_alerts_df = create_optimal_threshold_summary_df(:cases_after_alerts)
detected_outbreak_size_df = create_optimal_threshold_summary_df(
    :detected_outbreak_size
)
missed_outbreak_size_df = create_optimal_threshold_summary_df(
    :missed_outbreak_size
)

detectiondelays_wide_dfs = create_all_wide_optimal_threshold_summary_dfs(
    detectiondelays_df
)
cases_after_alerts_wide_df = create_all_wide_optimal_threshold_summary_dfs(
    cases_after_alerts_df
)
detected_outbreak_size_wide_dfs = create_all_wide_optimal_threshold_summary_dfs(
    detected_outbreak_size_df
)
