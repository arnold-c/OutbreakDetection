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

elisa_test_spec_vec = [
    # IndividualTestSpecification(0.98, 0.98, 0),
    IndividualTestSpecification(0.98, 0.98, 3),
    # IndividualTestSpecification(0.98, 0.98, 7),
    IndividualTestSpecification(0.98, 0.98, 14),
]

perfect_test_spec_vec = [
    IndividualTestSpecification(1.0, 1.0, 0),
    # IndividualTestSpecification(1.0, 1.0, 3),
    # IndividualTestSpecification(1.0, 1.0, 7),
    IndividualTestSpecification(1.0, 1.0, 14),
]

optimal_threshold_test_spec_vec = vcat(
    rdt_test_spec_vec, elisa_test_spec_vec, perfect_test_spec_vec
)

optimal_threshold_alertthreshold_vec = collect(1:1:15)

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

alert_method_vec = ["movingavg"]

#%%
ensemble_noise_specification =
    filter(ensemble_noise_specification_vec) do noise_spec
        noise_spec.noise_type == "poisson" ||
            noise_spec.correlation == "in-phase"
    end

ensemble_specification = ensemble_spec_vec[1]
alertmethod = alert_method_vec[1]

optimal_threshold_core_params = (
    alertthreshold_vec = optimal_threshold_alertthreshold_vec,
    ensemble_specification = ensemble_specification,
    outbreak_specification = ensemble_outbreak_specification,
    moving_avg_detection_lag = ensemble_moving_avg_detection_lag,
    percent_visit_clinic = ensemble_percent_visit_clinic,
    alertmethod = alertmethod,
)

clinical_hline = false

#%%
optimal_threshold_characteristics = collect_OptimalThresholdCharacteristics(
    ensemble_noise_specification,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    clinical_hline = clinical_hline,
);

#%%
function filter_optimal_threshold_characteristics_by_noise(
    optimal_threshold_characteristics,
    noise_specification,
)
    noise_optimal_solutions = filter(
        chars ->
            chars.noise_specification[1] == noise_specification,
        vec(optimal_threshold_characteristics),
    )
    @assert length(noise_optimal_solutions) == 1
    return noise_optimal_solutions[1]
end

#%%
poisson_noise_optimal_solutions = filter_optimal_threshold_characteristics_by_noise(
    optimal_threshold_characteristics,
    PoissonNoiseSpecification("poisson", 8.0);
);

#%%
dynamical_noise_optimal_solutions = filter_optimal_threshold_characteristics_by_noise(
    optimal_threshold_characteristics,
    DynamicalNoiseSpecification(
        "dynamical", 5.0, 7, 14, "in-phase", 0.15, 0.0492
    ),
);
