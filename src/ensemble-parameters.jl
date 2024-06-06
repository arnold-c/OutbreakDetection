using OutbreakDetectionUtils

ensemble_model_type = ("seasonal-infectivity-import", "tau-leaping")

nyears = 100
ensemble_time_specification = SimTimeParameters(;
    tmin = 0.0, tmax = 365.0 * nyears, tstep = 1.0
)

ensemble_state_specification = StateParameters(
    500_000,
    Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95),
)

ensemble_dynamics_specification = DynamicsParameters(
    ensemble_state_specification.init_states.N,
    27,
    0.2,
    SIGMA,
    GAMMA,
    16.0,
    0.8,
)

ensemble_nsims = 100

ensemble_specification = EnsembleSpecification(
    ensemble_model_type,
    ensemble_state_specification,
    ensemble_dynamics_specification,
    ensemble_time_specification,
    ensemble_nsims,
)

poisson_noise_mean_scaling_vec = [1.0, collect(2.0:2.0:8.0)...]

poisson_noise_spec_vec = create_combinations_vec(
    PoissonNoiseSpecification,
    (["poisson"], poisson_noise_mean_scaling_vec),
)

dynamical_noise_R0 = [5.0]
dynamical_noise_latent_period = [7]
dynamical_noise_duration_infection = [14]
dynamical_noise_correlation = ["in-phase"]#, "out-of-phase", "none"]
dynamical_noise_mean_scaling_vec = [0.15]
dynamical_noise_vaccination_coverage_vec = collect(0.0:0.2:0.8)
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

ensemble_noise_specification_vec = vcat(
    poisson_noise_spec_vec, dynamical_noise_spec_vec
)

ensemble_outbreak_specification = OutbreakSpecification(5, 30, 500)

ensemble_moving_avg_detection_lag = 7
ensemble_percent_visit_clinic = 0.6
ensemble_percent_clinic_tested_vec = collect(0.1:0.1:0.6)
