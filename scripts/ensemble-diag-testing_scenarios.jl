#%%
using DrWatson
@quickactivate "OutbreakDetection"

# include("../src/OutbreakDetection.jl")
# using OutbreakDetection

include("ensemble-sim_single-scenario_noise.jl")

#%%
model_types_vec = [("seasonal-infectivity-import", "tau-leaping")]
N_vec = [500_000]
nsims_vec = [1_000]
init_states_prop_dict = [Dict(:s_prop => 0.1, :e_prop => 0.01, :i_prop => 0.01, :r_prop => 0.88)]
tmin_vec = [0.0]
tstep_vec = [1.0]
tmax_vec = [365.0 * 100]

time_p_vec = vec(map(Iterators.product(tmin_vec, tstep_vec, tmax_vec)) do (tmin, tstep, tmax)
    SimTimeParameters(; tmin = tmin, tmax = tmax, tstep = tstep)
end)

beta_force_vec = [0.2]
# births_per_k_min = 5
# births_per_k_max = 20
# births_per_k_step = 5
# births_per_k_vec = collect(births_per_k_min:births_per_k_step:births_per_k_max)
births_per_k_vec = [5]
seed = 1234

ensemble_spec_vec = create_combinations_vec(
    EnsembleSpecification,
    (model_types_vec, N_vec, init_states_prop_dict, nsims_vec, births_per_k_vec, beta_force_vec, time_p_vec),
)

#%%
outbreak_threshold_vec = [5]
min_outbreak_dur_vec = [30]
min_outbreak_size_vec = [500]

outbreak_spec_vec = create_combinations_vec(
    OutbreakSpecification,
    (outbreak_threshold_vec, min_outbreak_dur_vec, min_outbreak_size_vec),
)

#%%
noise_spec_vec = [NoiseSpecification("static", ensemble_noise_arr)]

#%%
detectthreshold_vec = collect(5:5:20)
moveavglag_vec = [7]
perc_clinic_vec = [0.3]
perc_clinic_test_vec = [0.3]
testlag_vec = [3]

outbreak_detection_spec_vec = create_combinations_vec(
    OutbreakDetectionSpecification,
    (detectthreshold_vec, moveavglag_vec, perc_clinic_vec, perc_clinic_test_vec, testlag_vec),
)

#%%
testsens_vec = collect(0.8:0.1:1.0)
testspec_vec = collect(0.8:0.1:1.0)

test_spec_vec = create_combinations_vec(
    IndividualTestSpecification,
    (testsens_vec, testspec_vec),
)

#%%
base_scenarios_dict = @dict(
    ensemble_spec = ensemble_spec_vec,
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

