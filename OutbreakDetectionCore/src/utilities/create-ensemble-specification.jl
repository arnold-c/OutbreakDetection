export create_ensemble_specifications

"""
    create_ensemble_specifications(
        label,
		state_specification,
        time_specification,
        target_disease_dynamics_params::TargetDiseaseDynamicsParameters,
        common_disease_dynamics_params::CommonDiseaseDynamicsParameters,
        dynamical_noise_params::DynamicalNoiseParameters,
    )

Create a complete ensemble specification for disease simulation with emergent and null scenarios.

This function constructs an `EnsembleSpecification` by combining user-specified high-level
disease parameters with calculated epidemiological parameters. It creates both emergent
(outbreak) and null (baseline) dynamics specifications, where the null scenario maintains
the same disease dynamics but removes post-burnin vaccination coverage variation.

The function performs several key calculations:

  - Converts user-friendly disease parameters (R₀, latent period, infectious duration) into
    simulation-ready rates (β, σ, γ, μ)
  - Calculates the import rate (ε) based upon the formula presented on p210 of Keeling & Rohani
  - Determines minimum burnin vaccination coverage to achieve target effective reproduction number during the burnin period
  - Creates paired emergent and null dynamics specifications for comparative analysis

# Arguments

  - `label::String`: Disease name identifier for the ensemble
  - `time_specification`: Simulation time parameters including burnin duration and timestep
  - `state_specification`: Initial population state configuration (S, E, I, R compartments)
  - `target_disease_dynamics_params::TargetDiseaseDynamicsParameters`: High-level disease characteristics including R₀, disease periods, seasonality, and vaccination coverage ranges
  - `common_disease_dynamics_params::CommonDiseaseDynamicsParameters`: Shared parameters including birth rate, number of simulations, and burnin target Reff
  - `dynamical_noise_params::DynamicalNoiseParameters`: Specification for stochastic noise in disease parameters

# Returns

  - `EnsembleSpecification`: Complete specification ready for ensemble simulation, containing both emergent and null dynamics configurations

# Example

```julia
# Define disease parameters
target_params = TargetDiseaseDynamicsParameters(;
    R_0 = 2.5,
    latent_period_days = 5.0,
    infectious_duration_days = 7.0,
    beta_force = 0.1,
    seasonality = CosineSeasonality(),
    min_vaccination_coverage = 0.0,
    max_vaccination_coverage = 0.6,
)

common_params = CommonDiseaseDynamicsParameters(;
    births_per_k_pop = 14.0,
    nsims = 1000,
)

# Create ensemble specification
ensemble_spec = create_ensemble_specifications(
    "measles",
	state_spec,
    time_spec,
    target_params,
    common_params,
    noise_spec,
)
```

# See Also

  - [`EnsembleSpecification`](@ref): The returned ensemble configuration struct
  - [`TargetDiseaseDynamicsParameters`](@ref): High-level disease parameter specification
  - [`CommonDiseaseDynamicsParameters`](@ref): Shared simulation parameters
  - [`DynamicalNoiseParameters`](@ref): Stochastic noise specification
  - [`calculate_vaccination_rate_to_achieve_Reff`](@ref): Function for calculating vaccination coverage
"""
function create_ensemble_specifications(
        label,
        state_specification,
        time_specification,
        target_disease_dynamics_params::TargetDiseaseDynamicsParameters,
        common_disease_dynamics_params::CommonDiseaseDynamicsParameters,
        dynamical_noise_params::DynamicalNoiseParameters,
    )
    mu = calculate_mu(common_disease_dynamics_params.births_per_k_pop)
    gamma = calculate_gamma(target_disease_dynamics_params)
    sigma = calculate_sigma(target_disease_dynamics_params)

    beta_mean = calculate_beta(
        target_disease_dynamics_params.R_0,
        sigma,
        gamma,
        mu,
    )

    epsilon = calculate_import_rate(
        mu,
        target_disease_dynamics_params.R_0,
        state_specification.init_states.N,
    )

    dynamics_specification = DynamicsParameterSpecification(
        beta_mean,
        target_disease_dynamics_params.beta_force,
        target_disease_dynamics_params.seasonality,
        sigma,
        gamma,
        mu,
        common_disease_dynamics_params.births_per_k_pop,
        epsilon,
        target_disease_dynamics_params.R_0,
        target_disease_dynamics_params.min_vaccination_coverage,
        target_disease_dynamics_params.max_vaccination_coverage,
    )

    ensemble_specification = EnsembleSpecification(
        label,
        state_specification,
        time_specification,
        dynamics_specification,
        dynamical_noise_params,
        common_disease_dynamics_params.nsims,
    )

    return ensemble_specification
end
