export create_noise_dynamics_parameters

"""
    create_noise_dynamics_parameters(noise_spec, base_dynamics, population_N)

Create dynamics parameters for noise simulation.

This function creates a DynamicsParameters instance for noise simulations,
adjusting seasonality and transmission parameters based on the correlation
type specified in the noise specification.

# Arguments
- `noise_spec::DynamicalNoise`: Noise specification with vaccination coverage
- `base_dynamics::DynamicsParameterSpecification`: Base disease dynamics
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
noise_spec = DynamicalNoise(
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
base_dynamics = DynamicsParameterSpecification(target)

# Create noise dynamics
noise_dynamics = create_noise_dynamics_parameters(
    noise_spec,
    base_dynamics,
    500_000
)

# noise_dynamics has:
# - Same seasonality as base (in-phase)
# - R_0 = 5.0 (from noise_spec)
# - Different latent/infectious periods
# - vaccination_coverage = 0.65
```

# See Also
- [`DynamicalNoise`](@ref): Noise specification type
- [`DynamicsParameterSpecification`](@ref): Base dynamics type
- [`SeasonalityFunction`](@ref): Seasonality sum type
"""
function create_noise_dynamics_parameters(
        noise_spec::DynamicalNoise,
        base_dynamics::DynamicsParameterSpecification,
        population_N::Int64,
    )
    # Adjust beta_force based on correlation
    noise_beta_force = if noise_spec.correlation == "none"
        0.0
    else
        base_dynamics.beta_force
    end

    # Adjust seasonality based on correlation (using sum types)
    noise_seasonality = if noise_spec.correlation == "out-of-phase"
        # Flip seasonality for out-of-phase correlation
        base_variant = LightSumTypes.variant(base_dynamics.seasonality)
        if base_variant isa CosineSeasonality
            SeasonalityFunction(SineSeasonality())
        elseif base_variant isa SineSeasonality
            SeasonalityFunction(CosineSeasonality())
        else
            base_dynamics.seasonality
        end
    else
        base_dynamics.seasonality
    end

    # Calculate noise-specific parameters
    noise_gamma = 1.0 / noise_spec.duration_infection
    noise_sigma = 1.0 / noise_spec.latent_period
    noise_beta_mean = calculate_beta(
        noise_spec.R_0, noise_gamma, base_dynamics.mu, 1, population_N
    )
    noise_epsilon = calculate_import_rate(base_dynamics.mu, noise_spec.R_0, population_N)

    return DynamicsParameters(
        beta_mean = noise_beta_mean,
        beta_force = noise_beta_force,
        seasonality = noise_seasonality,
        sigma = noise_sigma,
        gamma = noise_gamma,
        mu = base_dynamics.mu,
        annual_births_per_k = base_dynamics.annual_births_per_k,
        epsilon = noise_epsilon,
        R_0 = noise_spec.R_0,
        vaccination_coverage = noise_spec.vaccination_coverage,
    )
end
