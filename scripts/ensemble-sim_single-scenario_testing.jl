#%%
using DrWatson
@quickactivate "OutbreakDetection"

include("ensemble-sim_single-scenario_noise.jl")

#%%
testlag = 3
moveavglag = 7
detectthreshold = 10
perc_clinic = 0.3
perc_clinic_test = 0.8
testsens = 0.9
testspec = 0.9

testing_arr = zeros(Int64, time_parameters.tlength, 8, size(inc_infec_arr, 3));
post_odds_arr = zeros(Float64, time_parameters.tlength, 2, size(inc_infec_arr, 3));

#%%
create_testing_arr!(
    testing_arr,
    inc_infec_arr,
    ensemble_noise_arr,
    post_odds_arr,
    OutbreakDetectionSpecification(detectthreshold, moveavglag, perc_clinic, perc_clinic_test, testlag),
    IndividualTestSpecification(testsens, testspec),
)

#%%
single_scenario_OT_Chars = calculate_OutbreakThresholdChars(testing_arr, inc_infec_arr)
