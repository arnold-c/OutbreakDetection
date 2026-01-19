export recreate_noise_dynamics_spec

"""
    create_noise_dynamics_parameters(noise_specification, dynamics_parameter_specification, population_N)

Create dynamics parameters for noise simulation.

This function creates a DynamicsParameters instance for noise simulations,
adjusting seasonality and transmission parameters based on the correlation
type specified in the noise specification.

# Arguments
- `noise_specification::DynamicalNoise`: Noise specification with vaccination coverage
- `dynamics_parameter_specification::DynamicsParameterSpecification`: Base disease dynamics
- `population_N::Int64`: Population size

# Returns
- `DynamicsParameters`: Dynamics parameters configured for noise simulation

# Correlation Types
- **"in-phase"**: Same seasonality as base dynamics, with seasonal forcing
- **"out-of-phase"**: Opposite seasonality (cosine ↔ sine), with seasonal forcing
- **"none"**: No seasonal forcing (beta_force = 0)

# Implementation Details
The function:
1. Adjusts beta_force based on correlation (0 for "none", base value otherwise)
2. Flips seasonality for "out-of-phase" (cosine ↔ sine)
3. Calculates noise-specific transmission parameters from noise R_0
4. Uses noise-specific latent period and infectious duration
5. Inherits demographic parameters (mu, births) from base dynamics

# Examples
```julia
# Create noise specification
noise_specification = DynamicalNoise(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 1.0,
    vaccination_coverage = 0.65
)

# Create base dynamics
target = TargetDiseaseDynamicsParameters(
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2
)
dynamics_parameter_specification = DynamicsParameterSpecification(target)

# Create noise dynamics
noise_dynamics = create_noise_dynamics_parameters(
    noise_specification,
    dynamics_parameter_specification,
    500_000
)

# noise_dynamics has:
# - Same seasonality as base (in-phase)
# - R_0 = 5.0 (from noise_specification)
# - Different latent/infectious periods
# - vaccination_coverage = 0.65
```

# See Also
- [`DynamicalNoise`](@ref): Noise specification type
- [`DynamicsParameterSpecification`](@ref): Base dynamics type
- [`SeasonalityFunction`](@ref): Seasonality sum type
"""
function recreate_noise_dynamics_spec(
        noise_specification::DynamicalNoise,
        ensemble_specification::EnsembleSpecification
    )
    UnPack.@unpack state_parameters,
        dynamics_parameter_specification = ensemble_specification

    N = state_parameters.init_states.N

    # Adjust beta_force based on correlation
    noise_beta_force = if noise_specification.correlation == "none"
        0.0
    else
        dynamics_parameter_specification.beta_force
    end

    # Adjust seasonality based on correlation (using sum types)
    noise_seasonality = if noise_specification.correlation == "out-of-phase"
        # Flip seasonality for out-of-phase correlation
        base_variant = LightSumTypes.variant(dynamics_parameter_specification.seasonality)
        if base_variant isa CosineSeasonality
            SeasonalityFunction(SineSeasonality())
        elseif base_variant isa SineSeasonality
            SeasonalityFunction(CosineSeasonality())
        else
            dynamics_parameter_specification.seasonality
        end
    else
        dynamics_parameter_specification.seasonality
    end

    # Calculate noise-specific parameters
    noise_gamma = calculate_gamma(dynamical_noise_params)
    noise_sigma = calculate_sigma(dynamical_noise_params)

    noise_beta_mean = calculate_beta(
        noise_specification.R_0,
        noise_sigma,
        noise_gamma,
        dynamics_parameter_specification.mu,
    )

    noise_epsilon = calculate_import_rate(
        dynamics_parameter_specification.mu,
        noise_specification.R_0,
        N
    )

    return DynamicsParameters(
        beta_mean = noise_beta_mean,
        beta_force = noise_beta_force,
        seasonality = noise_seasonality,
        sigma = noise_sigma,
        gamma = noise_gamma,
        mu = dynamics_parameter_specification.mu,
        annual_births_per_k = dynamics_parameter_specification.annual_births_per_k,
        epsilon = noise_epsilon,
        R_0 = noise_specification.R_0,
        vaccination_coverage = noise_specification.vaccination_coverage,
    )
end
