#%%
using DrWatson
@quickactivate "OutbreakDetection"

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(scriptsdir("ensemble-detection.jl"))
includet(scriptsdir("ensemble-noise-sim.jl"))
includet(funsdir("diag-testing-functions.jl"))

#%%
testlag = 3
moveavglag = 7
detectthreshold = 10
perc_clinic = 0.3
perc_clinic_test = 0.8
perc_tested = perc_clinic * perc_clinic_test
testsens = 0.9
testspec = 0.9

testing_arr = zeros(Int64, time_p.tlength, 8, size(inc_infec_arr, 3));
post_odds_arr = zeros(Float64, time_p.tlength, 2, size(inc_infec_arr, 3));

#%%
create_testing_arr!(
    testing_arr,
    inc_infec_arr,
    ensemble_noise_arr,
    post_odds_arr,
    perc_tested,
    testlag,
    testsens,
    testspec,
    detectthreshold,
    moveavglag,
)

#%%
OT_Chars = calculate_OutbreakThresholdChars(testing_arr, inc_infec_arr)

#%%
OT_Chars[1].crosstab
OT_Chars[1].outbreakbounds
OT_Chars[1].detectoutbreakbounds
OT_Chars[1].noutbreaks
OT_Chars[1].ndetectoutbreaks

# Note that an outbreak isn't detected continously!
testing_arr[80:100, :, 1]
