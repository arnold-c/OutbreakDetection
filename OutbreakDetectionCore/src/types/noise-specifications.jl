export NoiseType, PoissonNoiseType, DynamicalNoiseType,
    DynamicalNoiseSpecification, DynamicalNoise, PoissonNoise,
    NoiseSpecification, get_noise_description, get_noise_magnitude, getdirpath

# Abstract types for sum type variants
abstract type AbstractNoiseType end
struct PoissonNoiseType <: AbstractNoiseType end
struct DynamicalNoiseType <: AbstractNoiseType end

LightSumTypes.@sumtype NoiseType(PoissonNoiseType, DynamicalNoiseType) <: AbstractNoiseType

"""
    DynamicalNoiseSpecification

Specification for dynamical noise with optimization bounds.

This type defines the parameters for dynamical noise generation and the
bounds for vaccination coverage optimization. It does not include a specific
vaccination coverage value - that is determined by optimization and stored
in a DynamicalNoise instance.

# Fields
- `R_0::Float64`: Basic reproduction number for noise dynamics (must be > 0)
- `latent_period::Float64`: Latent period in days (must be > 0)
- `duration_infection::Float64`: Duration of infection in days (must be > 0)
- `correlation::String`: Correlation type ("in-phase", "out-of-phase", "none")
- `poisson_component::Float64`: Poisson noise component scaling (must be >= 0)
- `vaccination_bounds::Vector{Float64}`: Bounds for vaccination optimization [min, max]

# Validation
- R_0 must be positive
- Latent period must be positive
- Duration of infection must be positive
- Correlation must be one of: "in-phase", "out-of-phase", "none"
- Poisson component must be non-negative
- Vaccination bounds must have 2 elements with 0 <= min < max <= 1

# Examples
```julia
# In-phase correlation with Poisson component
spec = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 1.0,
    vaccination_bounds = [0.0, 0.8]
)

# Out-of-phase correlation without Poisson component
spec = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "out-of-phase",
    poisson_component = 0.0,
    vaccination_bounds = [0.0, 0.8]
)
```

# See Also
- [`DynamicalNoise`](@ref): Concrete instance with specific vaccination coverage
- [`PoissonNoise`](@ref): Simple Poisson noise alternative
"""
Base.@kwdef struct DynamicalNoiseSpecification
    R_0::Float64
    latent_period::Float64
    duration_infection::Float64
    correlation::String
    poisson_component::Float64
    vaccination_bounds::Vector{Float64} = [0.0, 1.0]

    function DynamicalNoiseSpecification(
            R_0, latent_period, duration_infection, correlation, poisson_component,
            vaccination_bounds
        )
        @assert R_0 > 0 "R_0 must be positive"
        @assert latent_period > 0 "Latent period must be positive"
        @assert duration_infection > 0 "Infectious duration must be positive"
        @assert correlation in ["in-phase", "out-of-phase", "none"] "Invalid correlation type"
        @assert poisson_component >= 0 "Poisson component must be non-negative"
        @assert length(vaccination_bounds) == 2 "Vaccination bounds must have 2 elements"
        @assert 0 <= vaccination_bounds[1] < vaccination_bounds[2] <= 1 "Invalid vaccination bounds"

        return new(
            R_0, latent_period, duration_infection, correlation, poisson_component,
            vaccination_bounds
        )
    end
end

"""
    DynamicalNoise

Concrete instance of dynamical noise with specific vaccination coverage.

This type represents a fully specified dynamical noise configuration,
including the vaccination coverage determined by optimization or manual
specification.

# Fields
- `R_0::Float64`: Basic reproduction number for noise dynamics
- `latent_period::Float64`: Latent period in days
- `duration_infection::Float64`: Duration of infection in days
- `correlation::String`: Correlation type ("in-phase", "out-of-phase", "none")
- `poisson_component::Float64`: Poisson noise component scaling
- `vaccination_coverage::Float64`: Specific vaccination coverage (0-1)

# Constructors
    DynamicalNoise(spec::DynamicalNoiseSpecification, vaccination_coverage::Float64)

Create a DynamicalNoise instance from a specification with a specific
vaccination coverage.

# Validation
- Vaccination coverage must be in [0, 1]

# Examples
```julia
spec = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 1.0,
    vaccination_bounds = [0.0, 0.8]
)

# Create instance with specific vaccination coverage
noise = DynamicalNoise(spec, 0.65)
```

# See Also
- [`DynamicalNoiseSpecification`](@ref): Specification with optimization bounds
- [`optimize_dynamic_noise_params`](@ref): Find optimal vaccination coverage
"""
Base.@kwdef struct DynamicalNoise
    R_0::Float64
    latent_period::Float64
    duration_infection::Float64
    correlation::String
    poisson_component::Float64
    vaccination_coverage::Float64
end

"""
    DynamicalNoise(spec::DynamicalNoiseSpecification, vaccination_coverage::Float64)

Create a DynamicalNoise instance from a specification.

# Arguments
- `spec::DynamicalNoiseSpecification`: Noise specification
- `vaccination_coverage::Float64`: Vaccination coverage (must be in [0, 1])

# Returns
- `DynamicalNoise`: Concrete noise instance

# Examples
```julia
spec = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 1.0
)

noise = DynamicalNoise(spec, 0.65)
```
"""
function DynamicalNoise(spec::DynamicalNoiseSpecification, vaccination_coverage::Float64)
    @assert 0.0 <= vaccination_coverage <= 1.0 "Vaccination coverage must be in [0, 1]"

    return DynamicalNoise(
        R_0 = spec.R_0,
        latent_period = spec.latent_period,
        duration_infection = spec.duration_infection,
        correlation = spec.correlation,
        poisson_component = spec.poisson_component,
        vaccination_coverage = vaccination_coverage,
    )
end

"""
    PoissonNoise

Simple Poisson-distributed noise specification.

This type represents noise generated purely from a Poisson distribution,
without any dynamical component. It is simpler than DynamicalNoise and
does not require optimization.

# Fields
- `noise_mean_scaling::Float64`: Scaling factor for noise mean (must be >= 0)

# Validation
- noise_mean_scaling must be non-negative

# Examples
```julia
# Low noise (1x scaling)
noise = PoissonNoise(noise_mean_scaling = 1.0)

# High noise (7x scaling)
noise = PoissonNoise(noise_mean_scaling = 7.0)
```

# See Also
- [`DynamicalNoise`](@ref): More complex dynamical noise alternative
"""
Base.@kwdef struct PoissonNoise
    noise_mean_scaling::Float64

    function PoissonNoise(noise_mean_scaling)
        @assert noise_mean_scaling >= 0 "Noise mean scaling must be non-negative"
        return new(noise_mean_scaling)
    end
end

# Sum type for polymorphic noise handling
abstract type AbstractNoiseSpecification end

LightSumTypes.@sumtype NoiseSpecification(PoissonNoise, DynamicalNoise) <:
AbstractNoiseSpecification

# Utility functions for backward compatibility

"""
    get_noise_description(noise::Union{PoissonNoise, DynamicalNoise})

Get a human-readable description of the noise type.

# Examples
```julia
poisson = PoissonNoise(noise_mean_scaling = 1.0)
get_noise_description(poisson)  # "Poisson"

dynamical = DynamicalNoise(spec, 0.65)
get_noise_description(dynamical)  # "Dynamical, in-phase"
```
"""
function get_noise_description(noise::PoissonNoise)
    return "Poisson"
end

function get_noise_description(noise::DynamicalNoise)
    return "Dynamical, $(noise.correlation)"
end

"""
    get_noise_magnitude(noise::Union{PoissonNoise, DynamicalNoise})

Get a description of the noise magnitude.

# Examples
```julia
poisson = PoissonNoise(noise_mean_scaling = 1.0)
get_noise_magnitude(poisson)  # "Poisson scaling: 1.0"

dynamical = DynamicalNoise(spec, 0.65)
get_noise_magnitude(dynamical)  # "Vaccination coverage: 0.65"
```
"""
function get_noise_magnitude(noise::PoissonNoise)
    return "Poisson scaling: $(noise.noise_mean_scaling)"
end

function get_noise_magnitude(noise::DynamicalNoise)
    return "Vaccination coverage: $(noise.vaccination_coverage)"
end

"""
    getdirpath(noise::Union{PoissonNoise, DynamicalNoise})

Generate directory path from noise properties.

# Examples
```julia
poisson = PoissonNoise(noise_mean_scaling = 1.0)
getdirpath(poisson)  # "noise_mean_scaling_1.0"

dynamical = DynamicalNoise(spec, 0.65)
getdirpath(dynamical)  # "R_0_5.0/latent_period_7.0/..."
```
"""
function getdirpath(noise::Union{PoissonNoise, DynamicalNoise})
    return reduce(
        joinpath, map(p -> "$(p)_$(getproperty(noise, p))", propertynames(noise))
    )
end

# Deprecated types for backward compatibility

"""
    PoissonNoiseSpecification

**DEPRECATED**: Use `PoissonNoise` instead.

This type is provided for backward compatibility but will be removed in a future version.
"""
struct PoissonNoiseSpecification
    noise_type::String
    noise_mean_scaling::Float64
end

function PoissonNoiseSpecification(noise_mean_scaling::Float64)
    @warn "PoissonNoiseSpecification is deprecated. Use PoissonNoise instead."
    return PoissonNoiseSpecification("Poisson", noise_mean_scaling)
end
