export DynamicsParameters,
    DynamicsParameterSpecification,
    TargetDiseaseDynamicsParameters,
    SeasonalityFunction,
    CosineSeasonality,
    SineSeasonality,
    CommonDiseaseDynamicsParameters

# Seasonality sum types
abstract type AbstractSeasonalityFunction end
struct CosineSeasonality <: AbstractSeasonalityFunction end
struct SineSeasonality <: AbstractSeasonalityFunction end

LightSumTypes.@sumtype SeasonalityFunction(
    CosineSeasonality, SineSeasonality
) <: AbstractSeasonalityFunction

"""
    TargetDiseaseDynamicsParameters

User-facing high-level specification of disease dynamics.

This is the primary interface for specifying disease parameters. Users provide
epidemiologically meaningful parameters, and the system calculates derived
quantities automatically.

# Fields

  - `R_0::Float64`: Basic reproduction number (must be > 0)
  - `latent_period_days::Float64`: Mean latent period in days (must be > 0)
  - `infectious_duration_days::Float64`: Mean infectious duration in days (must be > 0)
  - `beta_force::Float64`: Amplitude of seasonal forcing (must be in [0, 1])
  - `seasonality::SeasonalityFunction`: Seasonality function (default: CosineSeasonality)
  - `life_expectancy_years::Float64`: Life expectancy in years (default: 62.5)
  - `population_N::Int64`: Population size (default: 500,000)

# Examples

```julia
# Basic usage with defaults
target = TargetDiseaseDynamicsParameters(;
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2,
)

# With custom seasonality
target = TargetDiseaseDynamicsParameters(;
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2,
    seasonality = SeasonalityFunction(SineSeasonality()),
)
```

# See Also

  - [`DynamicsParameterSpecification`](@ref): Derived specification with calculated values
  - [`DynamicsParameters`](@ref): Concrete instance for simulation
"""
Base.@kwdef struct TargetDiseaseDynamicsParameters
    R_0::Float64
    latent_period::Dates.Day
    infectious_duration::Dates.Day
    beta_force::Float64
    seasonality::SeasonalityFunction = SeasonalityFunction(CosineSeasonality())
    min_vaccination_coverage::Float64 = 0.7
    max_vaccination_coverage::Float64 = 0.9

    function TargetDiseaseDynamicsParameters(
            R_0,
            latent_period,
            infectious_duration,
            beta_force,
            seasonality,
            min_vaccination_coverage,
            max_vaccination_coverage,
        )
        @assert R_0 > 0 "R_0 must be positive"
        @assert Dates.days(latent_period) > 0 "Latent period must be positive"
        @assert Dates.days(infectious_duration) > 0 "Infectious duration must be positive"
        @assert min_vaccination_coverage >= 0.0
        @assert min_vaccination_coverage <= max_vaccination_coverage <= 1.0

        return new(
            R_0,
            latent_period,
            infectious_duration,
            beta_force,
            seasonality,
            min_vaccination_coverage,
            max_vaccination_coverage,
        )
    end
end

Base.@kwdef struct CommonDiseaseDynamicsParameters
    births_per_k_pop::Float64
    nsims::Int64
end

"""
    DynamicsParameters

Concrete instance of disease dynamics parameters for simulation.

This is the final type used in actual simulations. It includes all parameters
from DynamicsParameterSpecification plus the vaccination coverage.

# Fields

  - `beta_mean::Float64`: Mean transmission rate
  - `beta_force::Float64`: Amplitude of seasonal forcing
  - `seasonality::SeasonalityFunction`: Seasonality function
  - `sigma::Float64`: Rate of progression from E to I (1/latent period)
  - `gamma::Float64`: Recovery rate (1/infectious period)
  - `mu::Float64`: Birth/death rate
  - `annual_births_per_k::Float64`: Annual births per 1000 population
  - `epsilon::Float64`: Import rate
  - `R_0::Float64`: Basic reproduction number
  - `vaccination_coverage::Float64`: Proportion vaccinated (0-1)

# Constructors

    DynamicsParameters(spec::DynamicsParameterSpecification; vaccination_coverage = 0.8)

Create concrete parameters from specification with specified vaccination coverage.

# Examples

```julia
# Full workflow
target = TargetDiseaseDynamicsParameters(;
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2,
)

spec = DynamicsParameterSpecification(target)
dynamics = DynamicsParameters(spec; vaccination_coverage = 0.8)
```

# See Also

  - [`TargetDiseaseDynamicsParameters`](@ref): User-facing specification
  - [`DynamicsParameterSpecification`](@ref): Derived specification
"""
Base.@kwdef struct DynamicsParameters
    beta_mean::Float64
    beta_force::Float64
    seasonality::SeasonalityFunction
    sigma::Float64
    gamma::Float64
    mu::Float64
    annual_births_per_k::Float64
    epsilon::Float64
    R_0::Float64
    vaccination_coverage::Float64
end

"""
    DynamicsParameterSpecification

Derived specification with calculated transmission and demographic parameters.

This intermediate type contains all calculated values needed for simulation,
derived from user-specified target parameters. It does not include
vaccination coverage, which is specified when creating a concrete
DynamicsParameters instance.

# Fields

  - `beta_mean::Float64`: Mean transmission rate
  - `beta_force::Float64`: Amplitude of seasonal forcing
  - `seasonality::SeasonalityFunction`: Seasonality function
  - `sigma::Float64`: Rate of progression from E to I (1/latent period)
  - `gamma::Float64`: Recovery rate (1/infectious period)
  - `mu::Float64`: Birth/death rate
  - `annual_births_per_k::Float64`: Annual births per 1000 population
  - `epsilon::Float64`: Import rate
  - `R_0::Float64`: Basic reproduction number
  - `population_N::Int64`: Population size
  - `min_vaccination_coverage::Float64`: The lower bound of the Uniform distribution to sample the vaccination coverage from
  - `max_vaccination_coverage::Float64`: The upper bound of the Uniform distribution to sample the vaccination coverage from

# Constructors

    DynamicsParameterSpecification(target::TargetDiseaseDynamicsParameters)

Create specification from target parameters with automatic calculation of
derived quantities.

# Examples

```julia
target = TargetDiseaseDynamicsParameters(;
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2,
)

spec = DynamicsParameterSpecification(target)
# spec.sigma ≈ 0.1 (1/10 days)
# spec.gamma ≈ 0.125 (1/8 days)
# spec.beta_mean calculated from R_0
```

# See Also

  - [`TargetDiseaseDynamicsParameters`](@ref): User-facing specification
  - [`DynamicsParameters`](@ref): Concrete instance for simulation
"""
Base.@kwdef struct DynamicsParameterSpecification
    beta_mean::Float64
    beta_force::Float64
    seasonality::SeasonalityFunction
    sigma::Float64
    gamma::Float64
    mu::Float64
    annual_births_per_k::Float64
    epsilon::Float64
    R_0::Float64
    min_vaccination_coverage::Float64
    max_vaccination_coverage::Float64
end

# Constructor: TargetDiseaseDynamicsParameters → DynamicsParameterSpecification
"""
    DynamicsParameterSpecification(target::TargetDiseaseDynamicsParameters)

Create a DynamicsParameterSpecification from target disease parameters.

Automatically calculates:

  - `sigma` from latent period
  - `gamma` from infectious duration
  - `mu` from life expectancy
  - `annual_births_per_k` from life expectancy
  - `beta_mean` from R_0 and other parameters
  - `epsilon` (import rate) from R_0 and population

# Validation

  - R_0 must be positive
  - Latent period must be positive
  - Infectious duration must be positive
  - Beta force must be in [0, 1]
  - Life expectancy must be positive
  - Population must be positive

# Examples

```julia
target = TargetDiseaseDynamicsParameters(;
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2,
)

spec = DynamicsParameterSpecification(target)
```
"""
function DynamicsParameterSpecification(
        state_specification::StateParameters,
        target_disease_dynamics_params::TargetDiseaseDynamicsParameters,
        common_disease_dynamics_params::CommonDiseaseDynamicsParameters,
    )
    @assert target_disease_dynamics_params.R_0 > 0 "R_0 must be positive"
    @assert Dates.days(target_disease_dynamics_params.latent_period) > 0 "Latent period must be positive"
    @assert Dates.days(target_disease_dynamics_params.infectious_duration) > 0 "Infectious duration must be positive"
    @assert 0 <= target_disease_dynamics_params.beta_force <= 1 "Beta force must be in [0, 1]"
    @assert state_specification.init_states.N > 0 "Population must be positive"

    sigma = calculate_sigma(target_disease_dynamics_params)
    gamma = calculate_gamma(target_disease_dynamics_params)
    mu = calculate_mu(common_disease_dynamics_params)

    beta_mean = calculate_beta(
        target_disease_dynamics_params.R_0,
        sigma,
        gamma,
        mu
    )

    epsilon = calculate_import_rate(
        mu,
        target_disease_dynamics_params.R_0,
        state_specification.init_states.N
    )

    return DynamicsParameterSpecification(;
        beta_mean = beta_mean,
        beta_force = target_disease_dynamics_params.beta_force,
        seasonality = target_disease_dynamics_params.seasonality,
        sigma = sigma,
        gamma = gamma,
        mu = mu,
        annual_births_per_k = common_disease_dynamics_params.births_per_k_pop,
        epsilon = epsilon,
        R_0 = target_disease_dynamics_params.R_0,
        min_vaccination_coverage = target_disease_dynamics_params.min_vaccination_coverage,
        max_vaccination_coverage = target_disease_dynamics_params.max_vaccination_coverage,
    )
end

# Constructor: DynamicsParameterSpecification → DynamicsParameters
"""
    DynamicsParameters(spec::DynamicsParameterSpecification; vaccination_coverage = 0.8)

Create concrete DynamicsParameters from specification with vaccination coverage.

# Arguments

  - `spec::DynamicsParameterSpecification`: Specification with calculated parameters
  - `vaccination_coverage::Float64`: Proportion vaccinated (default: 0.8, must be in [0, 1])

# Validation

  - Vaccination coverage must be in [0, 1]

# Examples

```julia
spec = DynamicsParameterSpecification(target)
dynamics = DynamicsParameters(spec)
```
"""
function DynamicsParameters(
        dynamic_parameter_specification::DynamicsParameterSpecification;
        seed = 1234,
    )
    Random.seed!(seed)

    vaccination_coverage =
    if dynamic_parameter_specification.min_vaccination_coverage ==
            dynamic_parameter_specification.max_vaccination_coverage
        dynamic_parameter_specification.min_vaccination_coverage
    else
        sample_vaccination_coverage(
            dynamic_parameter_specification.min_vaccination_coverage,
            dynamic_parameter_specification.max_vaccination_coverage,
        )
    end

    dynamics_parameters = DynamicsParameters(;
        beta_mean = dynamic_parameter_specification.beta_mean,
        beta_force = dynamic_parameter_specification.beta_force,
        seasonality = dynamic_parameter_specification.seasonality,
        sigma = dynamic_parameter_specification.sigma,
        gamma = dynamic_parameter_specification.gamma,
        mu = dynamic_parameter_specification.mu,
        annual_births_per_k = dynamic_parameter_specification.annual_births_per_k,
        epsilon = dynamic_parameter_specification.epsilon,
        R_0 = dynamic_parameter_specification.R_0,
        vaccination_coverage = vaccination_coverage,
    )
    return dynamics_parameters
end
