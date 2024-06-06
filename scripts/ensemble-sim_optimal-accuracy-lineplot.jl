#%%
using DrWatson
@quickactivate "OutbreakDetection"

using OutbreakDetectionUtils
using OutbreakDetection
using Revise

includet(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

if false
    include("../src/ensemble-parameters.jl")
end

#%%
optimal_threshold_test_spec_vec = [
    IndividualTestSpecification(0.85, 0.85, 0),
    # IndividualTestSpecification(0.9, 0.9, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
    # IndividualTestSpecification(1.0, 1.0, 3),
    # IndividualTestSpecification(1.0, 1.0, 7),
    # IndividualTestSpecification(1.0, 1.0, 14),
]

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
# for (ensemble_noise_specification, ensemble_specification, alertmethod) in
#     Iterators.product(
#     ensemble_noise_specification_vec[1], ensemble_spec_vec, alert_method_vec
# )
#     @info "Creating plots and tables for R0: $(ensemble_specification.dynamics_parameters.R_0), $(getdirpath(ensemble_noise_specification)), $(alertmethod)"
#     println("==============================================")

ensemble_noise_specification = ensemble_noise_specification_vec
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

# noise_specification_path = getdirpath(ensemble_noise_specification)
# noisespec_alertmethod_path = joinpath(noise_specification_path, alertmethod)
#
# if alertmethod != "dailythreshold"
#     noisespec_alertmethod_path = joinpath(
#         noisespec_alertmethod_path,
#         "moveavglag_$(ensemble_moving_avg_detection_lag)",
#     )
# end
#
# noisespec_alertmethod_filename = replace(
#     noisespec_alertmethod_path,
#     "/" => "_",
# )

basedirpath = joinpath(
    "R0_$(ensemble_specification.dynamics_parameters.R_0)"
    # noisespec_alertmethod_path,
)

baseplotdirpath = joinpath(
    plotsdir("ensemble/optimal-thresholds/lineplot"),
    basedirpath,
)

#%%
line_accuracy_plot(
    ensemble_noise_specification,
    ensemble_percent_clinic_tested_vec,
    optimal_threshold_test_spec_vec,
    optimal_threshold_core_params;
    plotdirpath = baseplotdirpath,
)
