export DynamicsParameters

"""
    DynamicsParameters

Parameters defining the disease dynamics and transmission.

# Fields
- `beta_mean::AbstractFloat`: Mean transmission rate
- `beta_force::AbstractFloat`: Amplitude of seasonal forcing
- `seasonality::Function`: Seasonality function (default: cos)
- `sigma::AbstractFloat`: Rate of progression from E to I (1/latent period)
- `gamma::AbstractFloat`: Recovery rate (1/infectious period)
- `mu::AbstractFloat`: Birth/death rate
- `annual_births_per_k::Union{Integer, AbstractFloat}`: Annual births per 1000
- `epsilon::AbstractFloat`: Import rate
- `R_0::AbstractFloat`: Basic reproduction number
- `vaccination_coverage::AbstractFloat`: Proportion vaccinated

# Constructors
Multiple constructors available for different use cases.
"""
struct DynamicsParameters{
        T1 <: AbstractFloat, T2 <: Union{<:Integer, T1}, T3 <: Function,
    }
    beta_mean::T1
    beta_force::T1
    seasonality::T3
    sigma::T1
    gamma::T1
    mu::T1
    annual_births_per_k::T2
    epsilon::T1
    R_0::T1
    vaccination_coverage::T1
end

function DynamicsParameters(
        sigma::Float64,
        gamma::Float64,
        R_0::Float64;
        vaccination_coverage::Float64 = 0.8,
    )
    return DynamicsParameters(
        BETA_MEAN,
        BETA_FORCE,
        cos,
        sigma,
        gamma,
        MU,
        ANNUAL_BIRTHS_PER_K,
        EPSILON,
        R_0,
        vaccination_coverage,
    )
end

function DynamicsParameters(
        N::Int64,
        annual_births_per_k::Int64,
        beta_force::Float64,
        sigma::Float64,
        gamma::Float64,
        R_0::Float64,
        vaccination_coverage::Float64;
        seasonality::Function = cos,
    )
    mu = calculate_mu(annual_births_per_k)
    beta_mean = calculate_beta(R_0, gamma, mu, 1, N)
    epsilon = calculate_import_rate(mu, R_0, N)

    return DynamicsParameters(
        beta_mean,
        beta_force,
        seasonality,
        sigma,
        gamma,
        mu,
        annual_births_per_k,
        epsilon,
        R_0,
        vaccination_coverage,
    )
end

function DynamicsParameters(
        N::Int64, annual_births_per_k::Int64, beta_force::Float64;
        vaccination_coverage::Float64 = 0.8,
    )
    mu = calculate_mu(annual_births_per_k)
    beta_mean = calculate_beta(R0, GAMMA, mu, 1, N)
    epsilon = calculate_import_rate(mu, R0, N)

    return DynamicsParameters(
        beta_mean,
        beta_force,
        cos,
        SIGMA,
        GAMMA,
        mu,
        annual_births_per_k,
        epsilon,
        R0,
        vaccination_coverage,
    )
end
