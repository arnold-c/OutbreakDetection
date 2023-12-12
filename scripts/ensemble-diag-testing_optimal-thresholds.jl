#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using FLoops
using NaNMath: NaNMath
using CSV: CSV
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
nyears::Float64 = 100.0
population_size::Int64 = 500_000

ensemble_specification = EnsembleSpecification(
    ("seasonal-infectivity-import", "tau-leaping"),
    StateParameters(
        population_size,
        Dict(
            :s_prop => 0.1,
            :e_prop => 0.0,
            :i_prop => 0.0,
            :r_prop => 0.9,
        ),
    ),
    DynamicsParameters(500_000, 10, 0.2; vaccination_coverage = 0.0),
    SimTimeParameters(;
        tmin = 0.0, tmax = 365.0 * nyears, tstep = 1.0
    ),
    100,
)
noise_specification = NoiseSpecification("poisson", 1.0)
outbreak_specification = OutbreakSpecification(5, 30, 500)

moving_avg_detection_lag = 7
percent_visit_clinic = 0.6
percent_clinic_tested_vec = collect(0.1:0.1:0.5)

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
            char = :unavoidable_cases,
            label = "Unavoidable Cases",
            color = (PERC_OUTBREAKS_MISSED_COLOR, 1.0),
            binwidth = 500,
        ),
        (
            char = :avoidable_cases,
            label = "Avoidable Cases",
            color = (PERC_OUTBREAKS_DETECTED_COLOR, 1.0),
            binwidth = 500,
        ),
    ],
)

#%%
cfr_df = CSV.read(
    datadir("CFR_2022.csv"),
    DataFrame; delim = ',',
    header = true,
    types = [String, Int64, Float64],
)
dropmissing!(cfr_df)

gha_cfr = only(cfr_df[cfr_df.country .== "GHA", :CFR])

population_df = CSV.read(
    datadir("input-populations.csv"),
    DataFrame; delim = ',',
    header = true
)

gha_2022_pop = only(population_df[population_df.ISO3_code .== "GHA", "2022"])
gha_2022_scale_population = gha_2022_pop / population_size

countries = [
    (;
        name = "Ghana",
        code = "GHA",
        cfr = gha_cfr,
        year = "2022",
        population_size = gha_2022_pop,
        scale_population = gha_2022_scale_population,
    ),
]

#%%
create_and_save_xlsx_optimal_threshold_summaries(optimal_thresholds_vec)

create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec, :detectiondelays
)

create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec,
    :unavoidable_cases;
    scale_annual = 1 / nyears,
    countries = countries,
)

create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec, :avoidable_cases;
    scale_annual = 1 / nyears,
    countries = countries,
)

create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec, :n_outbreak_cases;
    scale_annual = 1 / nyears,
    countries = countries,
)

create_and_save_xlsx_optimal_threshold_summaries(
    optimal_thresholds_vec, :n_tests
)
