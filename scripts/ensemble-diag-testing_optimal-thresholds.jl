#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using FLoops
using NaNMath: NaNMath
using DataFrames
using DataFramesMeta

using OutbreakDetection

includet(srcdir("makie-plotting-setup.jl"))

#%%
sensitivity_vec = collect(0.8:0.2:1.0)
specificity_vec = collect(0.8:0.2:1.0)
ind_test_spec_vec = Vector{IndividualTestSpecification}(
    undef,
    length(sensitivity_vec) + 1
)
for (i, sensitivity) in pairs(sensitivity_vec)
    ind_test_spec_vec[i] = IndividualTestSpecification(
        sensitivity, specificity_vec[i]
    )
end
ind_test_spec_vec[end] = IndividualTestSpecification(1.0, 0.0)

detectthreshold_vec = collect(4:1:15)

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
test_result_lag = 0
percent_visit_clinic = 0.6
percent_clinic_tested_vec = collect(0.2:0.2:1.0)

threshold_comparison_params = (
    detectthreshold_vec = detectthreshold_vec,
    ensemble_specification = ensemble_specification,
    noise_specification = noise_specification,
    outbreak_specification = outbreak_specification,
    moving_avg_detection_lag = moving_avg_detection_lag,
    test_result_lag = test_result_lag,
    percent_visit_clinic = percent_visit_clinic,
)

#%%
optimal_thresholds_vec = calculate_OptimalThresholdCharacteristics(
    percent_clinic_tested_vec,
    ind_test_spec_vec,
    threshold_comparison_params
)

#%%
testspecs = StructArray(
    getfield.(
        optimal_thresholds_vec.scenario_specification,
        :individual_test_specification,
    ),
)

clinictesting =
    getfield.(
        getfield.(
            optimal_thresholds_vec.scenario_specification,
            :outbreak_detection_specification,
        ),
        :percent_clinic_tested,
    )
optimal_thresholds_df = DataFrame(;
    clinictesting,
    sensitivity = testspecs.sensitivity,
    specificity = testspecs.specificity,
    detection_threshold = optimal_thresholds_vec.detection_threshold,
    accuracy = optimal_thresholds_vec.accuracy,
)

@chain optimal_thresholds_df begin
    @orderby :specificity
    unstack(
        _, [:sensitivity, :specificity], :clinictesting, :detection_threshold
    )
end

@chain optimal_thresholds_df begin
    @orderby :specificity
    unstack(
        _, [:sensitivity, :specificity], :clinictesting, :accuracy
    )
end
