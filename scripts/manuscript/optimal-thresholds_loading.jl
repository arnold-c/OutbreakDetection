#%%
using DrWatson
@quickactivate "OutbreakDetection"

using OutbreakDetectionUtils:
    IndividualTestSpecification, create_combinations_vec,
    DynamicalNoiseSpecification, PoissonNoiseSpecification
using OutbreakDetection: collect_OptimalThresholdCharacteristics

include(srcdir("ensemble-parameters.jl"))

if false
    include("../src/ensemble-parameters.jl")
end

#%%
rdt_test_spec_vec = [
    IndividualTestSpecification(0.85, 0.85, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
]

perfect_test_spec_vec = [
    IndividualTestSpecification(1.0, 1.0, 0),
    IndividualTestSpecification(1.0, 1.0, 14),
]

test_spec_vec = vcat(
    rdt_test_spec_vec,
    perfect_test_spec_vec,
)

#%%
R_0_vec = [16.0]

ensemble_dynamics_spec_vec = create_combinations_vec(
    DynamicsParameters,
    (
        [ensemble_state_specification.init_states.N],
        [27],
        [0.2],
        [SIGMA],
        [GAMMA],
        R_0_vec,
        [0.8],
    ),
)

ensemble_spec_vec = create_combinations_vec(
    EnsembleSpecification,
    (
        [ensemble_model_type],
        [ensemble_state_specification],
        ensemble_dynamics_spec_vec,
        [ensemble_time_specification],
        [ensemble_nsims],
    ),
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
noise_spec_vec =
    filter(ensemble_noise_specification_vec) do noise_spec
    noise_spec.noise_type == "poisson" ||
        noise_spec.correlation == "in-phase"
end

ensemble_specification = ensemble_spec_vec[1]
alertmethod = alert_method_vec[1]

clinical_hline = false

accuracy_functions = [arithmetic_mean, calculate_f_beta_score]

#%%
all_optim_df = OutbreakDetectionUtils.run_scenario_optimizations(
    ensemble_spec_vec,
    outbreak_spec_vec,
    noise_spec_vec,
    outbreak_detection_spec_vec,
    test_spec_vec,
    MSO,
    accuracy_functions;
    force = false,
    save_df = true,
    return_df = true,
    filter_df_results = true,
)

optim_df = DataFrames.subset(
    all_optim_df,
    :accuracy_function => f -> f .== arithmetic_mean,
)

f1_optim_df = DataFrames.subset(
    all_optim_df,
    :accuracy_function => f -> f .== calculate_f_beta_score,
)

#%%
poisson_noise_optimal_solutions = DataFrames.filter(
    :noise_spec => n -> n == PoissonNoiseSpecification("poisson", 8.0),
    optim_df,
);

dynamical_noise_optimal_solutions = DataFrames.filter(
    :noise_spec =>
        n ->
    n == DynamicalNoiseSpecification(
        "dynamical", 5.0, 7, 14, "in-phase", 0.15, 0.0492
    ),
    optim_df,
);

#%%
optimal_threshold_characteristics = reshape_optim_df_to_matrix(optim_df);

f1_optimal_threshold_characteristics = reshape_optim_df_to_matrix(f1_optim_df);
