#%%
using DrWatson
@quickactivate "OutbreakDetection"

using OutbreakDetectionUtils
using OutbreakDetection
using FLoops
using JLD2

include(srcdir("makie-plotting-setup.jl"))

include(srcdir("ensemble-parameters.jl"))

if false
    include("../src/ensemble-parameters.jl")
end

#%%
model_types_vec = [("seasonal-infectivity-import", "tau-leaping")]

#%%
N_vec = [500_000]
nsims_vec = [100]
init_states_prop_dict = [
    Dict(:s_prop => 0.05, :e_prop => 0.00, :i_prop => 0.00, :r_prop => 0.95)
]

ensemble_state_p_vec = create_combinations_vec(
    StateParameters, (N_vec, init_states_prop_dict)
)

#%%
tmin_vec = [0.0]
tstep_vec = [1.0]
tmax_vec = [365.0 * 100]

time_p_vec = vec(
    map(
        Iterators.product(tmin_vec, tstep_vec, tmax_vec)
    ) do (tmin, tstep, tmax)
        SimTimeParameters(; tmin = tmin, tmax = tmax, tstep = tstep)
    end,
)

#%%
beta_force_vec = [0.2]
annual_births_per_k_vec = [27]
seed = 1234

#%%
latent_per_days_vec = [LATENT_PER_DAYS]
dur_inf_days_vec = [DUR_INF_DAYS]
R_0_vec = [16.0]
sigma_vec = 1 ./ latent_per_days_vec
gamma_vec = 1 ./ dur_inf_days_vec
vaccination_coverage_vec = [0.8]

#%%
ensemble_spec_vec = create_ensemble_spec_combinations(
    beta_force_vec,
    [cos],
    sigma_vec,
    gamma_vec,
    annual_births_per_k_vec,
    R_0_vec,
    vaccination_coverage_vec,
    N_vec,
    init_states_prop_dict,
    model_types_vec,
    time_p_vec,
    nsims_vec,
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
# TODO: How to make the mean of the noise a proportion of the mean of incidence
# when incidence depends on the dynamics parameters?
# Could pass variables to ensemble function and calculate each simulations and
# scenario's noise mean, but that would break implementation using NoiseSpecification
# struct currently
poisson_noise_mean_scaling_vec = [
	1.0,
	collect(2.0:2.0:8.0)...
]

poisson_noise_spec_vec = create_combinations_vec(
    PoissonNoiseSpecification,
    (["poisson"], poisson_noise_mean_scaling_vec),
)

dynamical_noise_R0 = [5.0]
dynamical_noise_latent_period = [7]
dynamical_noise_duration_infection = [14]
dynamical_noise_correlation = [
    "in-phase"
    # "out-of-phase",
    # "none",
]
dynamical_noise_mean_scaling_vec = [0.15]
dynamical_noise_vaccination_coverage_vec = [
    0.8538, 0.7383, 0.5088, 0.2789, 0.0492
]
dynamical_noise_spec_vec = create_combinations_vec(
    DynamicalNoiseSpecification,
    (
        ["dynamical"],
        dynamical_noise_R0,
        dynamical_noise_latent_period,
        dynamical_noise_duration_infection,
        dynamical_noise_correlation,
        dynamical_noise_mean_scaling_vec,
        dynamical_noise_vaccination_coverage_vec,
    ),
)

noise_spec_vec = vcat(
	poisson_noise_spec_vec,
	dynamical_noise_spec_vec
)

#%%
# This will get updated in optimization so just use placeholder
alertthreshold_vec = 5.0
moveavglag_vec = [7]
perc_clinic_vec = [1.0]
perc_clinic_test_vec = collect(0.1:0.1:1.0)
alert_method_vec = ["movingavg"]

outbreak_detection_spec_vec = create_combinations_vec(
    OutbreakDetectionSpecification,
    (
        alertthreshold_vec,
        moveavglag_vec,
        perc_clinic_vec,
        perc_clinic_test_vec,
        alert_method_vec,
    ),
)

#%%
test_spec_vec = [
    # PROTOTYPE_RDT_TEST_SPECS...,
    IndividualTestSpecification(0.85, 0.85, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    # CLINICAL_TEST_SPECS...,
    # IndividualTestSpecification(0.98, 0.98, 0),
    # IndividualTestSpecification(0.98, 0.98, 3),
    # IndividualTestSpecification(0.98, 0.98, 7),
    # IndividualTestSpecification(0.98, 0.98, 14),
    IndividualTestSpecification(1.0, 1.0, 0),
    # IndividualTestSpecification(1.0, 1.0, 3),
    # IndividualTestSpecification(1.0, 1.0, 7),
    IndividualTestSpecification(1.0, 1.0, 14),
]

#%%
optim_df = OutbreakDetectionUtils.run_scenario_optimizations(
    ensemble_spec_vec,
    outbreak_spec_vec,
    noise_spec_vec,
    outbreak_detection_spec_vec,
    test_spec_vec,
    MSO;
    executor = SequentialEx(),
)

#%%
most_recent_optim_df = get_most_recent_optimization_filepath(
	"alert-threshold-optimization.jld2",
    outdir("ensemble", "scenario-optimizations"),
) |>
	fp -> Try.unwrap(fp) |>
	fp -> JLD2.load(fp)

most_recent_optim_df["threshold_optim_df"]

# #%%
# missing_optimizations = check_missing_scenario_optimizations(
# 	optim_df,
#     ensemble_spec_vec,
#     outbreak_spec_vec,
#     noise_spec_vec,
#     outbreak_detection_spec_vec,
#     test_spec_vec,
#     MSO;
# )
#
# #%%
# run_missing_scenario_optimizations!(
# 	optim_df,
# 	missing_optimizations
# )
#
# #%%
# optim_df
#
# #%%
# @tagsave(outdir("2025-03-13_11:00:00_optimization-df.jld2"), Dict("optim_df" => optim_df))
#
#
# #%%
# base_param_dict = @dict(
# 	ensemble_spec = ensemble_spec_vec[1],
# 	outbreak_spec = outbreak_spec_vec[1],
# 		seed = 1234,
# )
#
# ensemble_inc_arr, thresholds_vec = setup_optimization(
# 	base_param_dict
# )
#
# noise_array, noise_means = create_noise_arr(
# 	noise_spec_vec[1],
# 	ensemble_inc_arr;
# 	ensemble_specification = ensemble_spec_vec[1],
# 	seed = seed,
# )
#
# obj_inputs = (;
# 	ensemble_inc_arr,
# 	noise_array,
# 	outbreak_detection_specification = outbreak_detection_spec_vec[1],
# 	individual_test_specification = test_spec_vec[2],
# 	thresholds_vec,
# )
#
# objective_function_closure =
# 	x -> objective_function(x, obj_inputs)
#
# #%%
# @elapsed optim_minimizer, optim_minimum = optimization_wrapper(
# 	objective_function_closure,
# 	QD;
# )
# @show optim_minimizer, optim_minimum
#
# #%%
# @elapsed optim_minimizer, optim_minimum = optimization_wrapper(
# 	objective_function_closure,
# 	MSO;
# )
#
# @show optim_minimizer, optim_minimum
