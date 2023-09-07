#%%
using DrWatson
@quickactivate "OutbreakDetection"

# include("../src/OutbreakDetection.jl")
# using .OutbreakDetection

include("ensemble-sim_single-scenario_noise.jl")

#%%
outbreak_threshold_vec = [5]
min_outbreak_dur_vec = [30]
min_outbreak_size_vec = [500]

outbreak_spec_combinations = Iterators.product(
    outbreak_threshold_vec,
    min_outbreak_dur_vec,
    min_outbreak_size_vec
)

outbreak_spec_vec = vec(
    map(
        combination -> OutbreakSpecification(combination...),
        outbreak_spec_combinations,
    ),
)

#%%
noise_spec_vec = [NoiseSpecification("static", ensemble_noise_arr)]

#%%
detectthreshold_vec = collect(5:5:20)
moveavglag_vec = [7]
perc_clinic_vec = collect(0.3)
perc_clinic_test_vec = collect(0.3)
testlag_vec = collect(3)

outbreak_detection_spec_combinations = Iterators.product(
    detectthreshold_vec,
    moveavglag_vec,
    perc_clinic_vec,
    perc_clinic_test_vec,
    testlag_vec,
)

outbreak_detection_spec_vec = vec(
    map(
        combination -> OutbreakDetectionSpecification(combination...),
        outbreak_detection_spec_combinations,
    ),
)

#%%
testsens_vec = collect(0.8:0.1:1.0)
testspec_vec = collect(0.8:0.1:1.0)

test_spec_combinations = Iterators.product(
    testsens_vec,
    testspec_vec,
)

test_spec_vec = vec(
    map(
        combination -> IndividualTestSpecification(combination...),
        test_spec_combinations,
    ),
)

#%%
scenarios_dict = dict_list(
    @dict(
        outbreak_spec = outbreak_spec_vec,
        noise_spec = noise_spec_vec,
        outbreak_detect_spec = outbreak_detection_spec_vec,
        ind_test_spec = test_spec_vec,
    )
)

#%%

