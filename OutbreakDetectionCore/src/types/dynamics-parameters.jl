export TargetDiseaseDynamicsParameters, DynamicsParameterSpecification,
    DynamicsParameters, SeasonalityFunction, CosineSeasonality, SineSeasonality

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
target = TargetDiseaseDynamicsParameters(
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2
)

# With custom seasonality
target = TargetDiseaseDynamicsParameters(
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2,
    seasonality = SeasonalityFunction(SineSeasonality())
)
```

# See Also
- [`DynamicsParameterSpecification`](@ref): Derived specification with calculated values
- [`DynamicsParameters`](@ref): Concrete instance for simulation
"""
Base.@kwdef struct TargetDiseaseDynamicsParameters
    R_0::Float64
    latent_period_days::Float64
    infectious_duration_days::Float64
    beta_force::Float64
    seasonality::SeasonalityFunction = SeasonalityFunction(CosineSeasonality())
    life_expectancy_years::Float64 = 62.5
    population_N::Int64 = 500_000
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

# Constructors
    DynamicsParameterSpecification(target::TargetDiseaseDynamicsParameters)

Create specification from target parameters with automatic calculation of
derived quantities.

# Examples
```julia
target = TargetDiseaseDynamicsParameters(
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2
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
    population_N::Int64
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
target = TargetDiseaseDynamicsParameters(
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2
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
target = TargetDiseaseDynamicsParameters(
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2
)

spec = DynamicsParameterSpecification(target)
```
"""
function DynamicsParameterSpecification(target::TargetDiseaseDynamicsParameters)
    @assert target.R_0 > 0 "R_0 must be positive"
    @assert target.latent_period_days > 0 "Latent period must be positive"
    @assert target.infectious_duration_days > 0 "Infectious duration must be positive"
    @assert 0 <= target.beta_force <= 1 "Beta force must be in [0, 1]"
    @assert target.life_expectancy_years > 0 "Life expectancy must be positive"
    @assert target.population_N > 0 "Population must be positive"

    sigma = 1.0 / target.latent_period_days
    gamma = 1.0 / target.infectious_duration_days
    mu = 1.0 / (target.life_expectancy_years * 365.0)
    annual_births_per_k = 1000.0 / target.life_expectancy_years
    beta_mean = calculate_beta(target.R_0, gamma, mu, 1, target.population_N)
    epsilon = calculate_import_rate(mu, target.R_0, target.population_N)

    return DynamicsParameterSpecification(
        beta_mean = beta_mean,
        beta_force = target.beta_force,
        seasonality = target.seasonality,
        sigma = sigma,
        gamma = gamma,
        mu = mu,
        annual_births_per_k = annual_births_per_k,
        epsilon = epsilon,
        R_0 = target.R_0,
        population_N = target.population_N,
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
dynamics = DynamicsParameters(spec; vaccination_coverage = 0.8)
```
"""
function DynamicsParameters(
        spec::DynamicsParameterSpecification; vaccination_coverage::Float64 = 0.8
    )
    @assert 0.0 <= vaccination_coverage <= 1.0 "Vaccination coverage must be in [0, 1]"

    return DynamicsParameters(
        beta_mean = spec.beta_mean,
        beta_force = spec.beta_force,
        seasonality = spec.seasonality,
        sigma = spec.sigma,
        gamma = spec.gamma,
        mu = spec.mu,
        annual_births_per_k = spec.annual_births_per_k,
        epsilon = spec.epsilon,
        R_0 = spec.R_0,
        vaccination_coverage = vaccination_coverage,
    )
end

# Deprecated constructors for backward compatibility
"""
    DynamicsParameters(sigma, gamma, R_0; vaccination_coverage = 0.8)

**DEPRECATED**: Use TargetDiseaseDynamicsParameters → DynamicsParameterSpecification → DynamicsParameters

This constructor is provided for backward compatibility but will be removed in a future version.
"""
function DynamicsParameters(
        sigma::Float64, gamma::Float64, R_0::Float64; vaccination_coverage::Float64 = 0.8
    )
    @warn "DynamicsParameters(sigma, gamma, R_0) is deprecated. " *
        "Use TargetDiseaseDynamicsParameters → DynamicsParameterSpecification → DynamicsParameters"

    # Provide default values for migration
    target = TargetDiseaseDynamicsParameters(
        R_0 = R_0,
        latent_period_days = 1.0 / sigma,
        infectious_duration_days = 1.0 / gamma,
        beta_force = 0.2,  # default
    )
    spec = DynamicsParameterSpecification(target)
    return DynamicsParameters(spec; vaccination_coverage = vaccination_coverage)
end

"""
    DynamicsParameters(N, annual_births_per_k, beta_force, sigma, gamma, R_0, vaccination_coverage; seasonality = cos)

**DEPRECATED**: Use TargetDiseaseDynamicsParameters → DynamicsParameterSpecification → DynamicsParameters

This constructor is provided for backward compatibility but will be removed in a future version.
"""
function DynamicsParameters(
        N::Int64,
        annual_births_per_k::Union{Int64, Float64},
        beta_force::Float64,
        sigma::Float64,
        gamma::Float64,
        R_0::Float64,
        vaccination_coverage::Float64;
        seasonality::Union{Function, SeasonalityFunction} = cos,
    )
    @warn "DynamicsParameters(N, annual_births_per_k, ...) is deprecated. " *
        "Use TargetDiseaseDynamicsParameters → DynamicsParameterSpecification → DynamicsParameters"

    # Convert seasonality function to sum type if needed
    seasonality_sumtype = if seasonality isa Function
        if seasonality === cos
            SeasonalityFunction(CosineSeasonality())
        elseif seasonality === sin
            SeasonalityFunction(SineSeasonality())
        else
            @warn "Unknown seasonality function, defaulting to cosine"
            SeasonalityFunction(CosineSeasonality())
        end
    else
        seasonality
    end

    mu = calculate_mu(annual_births_per_k)
    beta_mean = calculate_beta(R_0, gamma, mu, 1, N)
    epsilon = calculate_import_rate(mu, R_0, N)

    return DynamicsParameters(
        beta_mean = beta_mean,
        beta_force = beta_force,
        seasonality = seasonality_sumtype,
        sigma = sigma,
        gamma = gamma,
        mu = mu,
        annual_births_per_k = annual_births_per_k,
        epsilon = epsilon,
        R_0 = R_0,
        vaccination_coverage = vaccination_coverage,
    )
end

"""
    DynamicsParameters(N, annual_births_per_k, beta_force; vaccination_coverage = 0.8)

**DEPRECATED**: Use TargetDiseaseDynamicsParameters → DynamicsParameterSpecification → DynamicsParameters

This constructor is provided for backward compatibility but will be removed in a future version.
"""
function DynamicsParameters(
        N::Int64,
        annual_births_per_k::Union{Int64, Float64},
        beta_force::Float64;
        vaccination_coverage::Float64 = 0.8,
    )
    @warn "DynamicsParameters(N, annual_births_per_k, beta_force) is deprecated. " *
        "Use TargetDiseaseDynamicsParameters → DynamicsParameterSpecification → DynamicsParameters"

    # Use constants from the old system
    # Note: This requires constants to still be available during transition
    target = TargetDiseaseDynamicsParameters(
        R_0 = 16.0,  # Default R0 from constants
        latent_period_days = 10.0,  # Default from constants
        infectious_duration_days = 8.0,  # Default from constants
        beta_force = beta_force,
        life_expectancy_years = 1000.0 / annual_births_per_k,
        population_N = N,
    )

    spec = DynamicsParameterSpecification(target)
    return DynamicsParameters(spec; vaccination_coverage = vaccination_coverage)
end
