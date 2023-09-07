#%%
using DrWatson
@quickactivate "OutbreakDetection"

# include("../src/OutbreakDetection.jl")
# using OutbreakDetection

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
perc_clinic_vec = [0.3]
perc_clinic_test_vec = [0.3]
testlag_vec = [3]

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
base_scenarios_dict = @dict(
    outbreak_spec = outbreak_spec_vec,
    noise_spec = noise_spec_vec,
    outbreak_detect_spec = outbreak_detection_spec_vec,
    ind_test_spec = test_spec_vec,
)

#%%
N_vec = convert.(Int64, [5e5])
nsims_vec = [1_000]
init_states_prop_dict = [
    Dict(
        :s_prop => 0.1,
        :e_prop => 0.01,
        :i_prop => 0.01,
        :r_prop => 0.88,
    )
]
tstep_vec = [1.0]
tmax_vec = [365.0 * 100]
beta_force_vec = [0.2]
births_per_k_vec = [10]

ensemble_time_p = SimTimeParameters(;
    tmin = 0.0, tmax = 365.0 * 100.0, tstep = 1.0
)

base_param_dict = @dict(
    N = N_vec,
    init_states_prop = init_states_prop_dict,
    time_p = ensemble_time_p,
    nsims = nsims_vec,
    beta_force = beta_force_vec,
    births_per_k = births_per_k_vec,
    seed = seed,
)

#%%
scenarios_dict = dict_list(merge(base_param_dict, base_scenarios_dict))

#%%
run_OutbreakThresholdChars_creation(scenarios_dict; progress = true)

#%%
scenario = ScenarioSpecification(
    ensemble_spec,
    outbreak_spec,
    NoiseSpecification("static", ensemble_noise_arr),
    OutbreakDetectionSpecification(10, 7, 0.3, 0.3, 3),
    IndividualTestSpecification(0.8, 0.8),
)

scenario_file = get_scenario_file(scenario)

@unpack OT_chars = scenario_file
