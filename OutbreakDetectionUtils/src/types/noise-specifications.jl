export NoiseSpecification, PoissonNoiseSpecification, DynamicalNoiseSpecification,
    get_noise_description, get_noise_magnitude, getdirpath

"""
    NoiseSpecification

Abstract type for noise specifications.
"""
abstract type NoiseSpecification end

"""
    PoissonNoiseSpecification

Specification for Poisson-distributed noise.

# Fields
- `noise_type::AbstractString`: Type of noise ("Poisson")
- `noise_mean_scaling::AbstractFloat`: Scaling factor for noise mean
"""
struct PoissonNoiseSpecification{
        T1 <: AbstractString, T2 <: AbstractFloat,
    } <: NoiseSpecification
    noise_type::T1
    noise_mean_scaling::T2
end

function PoissonNoiseSpecification(noise_mean_scaling)
    return PoissonNoiseSpecification("Poisson", noise_mean_scaling)
end

"""
    DynamicalNoiseSpecification

Specification for dynamically-generated noise.

# Fields
- `noise_type::AbstractString`: Type of noise
- `R_0::AbstractFloat`: Basic reproduction number for noise dynamics
- `latent_period::Integer`: Latent period in days
- `duration_infection::Integer`: Duration of infection in days
- `correlation::AbstractString`: Correlation type
- `noise_mean_scaling::AbstractFloat`: Scaling factor for noise mean
- `vaccination_coverage::AbstractFloat`: Vaccination coverage
"""
struct DynamicalNoiseSpecification{
        T1 <: AbstractString, T2 <: AbstractFloat, T3 <: Integer,
    } <: NoiseSpecification
    noise_type::T1
    R_0::T2
    latent_period::T3
    duration_infection::T3
    correlation::T1
    noise_mean_scaling::T2
    vaccination_coverage::T2
end

function DynamicalNoiseSpecification()
end

"""
    get_noise_description(noise_specification::NoiseSpecification)

Get a human-readable description of the noise specification.
"""
function get_noise_description(
        noise_specification::T
    ) where {T <: NoiseSpecification}
    return noise_specification.noise_type
end

function get_noise_description(noise_specification::DynamicalNoiseSpecification)
    return string(
        noise_specification.noise_type, ", ", noise_specification.correlation
    )
end

"""
    get_noise_magnitude(noise_specification::NoiseSpecification)

Get a description of the noise magnitude.
"""
function get_noise_magnitude(
        noise_specification::T
    ) where {T <: NoiseSpecification}
    return string("Poisson scaling: ", noise_specification.noise_mean_scaling)
end

function get_noise_magnitude(
        noise_specification::DynamicalNoiseSpecification
    )
    return string("Rubella vax: ", noise_specification.vaccination_coverage)
end

"""
    getdirpath(spec::NoiseSpecification)

Generate directory path from noise specification properties.
"""
function getdirpath(spec::NoiseSpecification)
    return reduce(
        joinpath,
        map(
            p -> "$(p)_$(getproperty(spec, p))",
            propertynames(spec),
        ),
    )
end
