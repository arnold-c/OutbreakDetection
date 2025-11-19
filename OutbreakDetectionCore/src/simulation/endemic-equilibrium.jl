export calculate_endemic_equilibrium_proportions

"""
    calculate_endemic_equilibrium_proportions(dynamics_params, vaccination_coverage)

Calculate endemic equilibrium proportions for SEIR model with vaccination.

This function computes the equilibrium proportions of susceptible, exposed,
infected, and recovered individuals for an SEIR model with vaccination at birth.
The equilibrium exists only when the effective reproduction number R_eff > 1.

# Arguments
- `dynamics_params::Union{DynamicsParameters, DynamicsParameterSpecification}`: Disease dynamics
- `vaccination_coverage::Float64`: Vaccination coverage (0-1)

# Returns
- `Try.Ok((s_prop, e_prop, i_prop, r_prop))` on success
- `Try.Err(message)` if equilibrium doesn't exist (R_eff ≤ 1) or parameters invalid

# Mathematical Background
For an SEIR model with vaccination coverage ρ:
- R_eff = R₀(1 - ρ)
- s* = 1/R₀
- i* = μ(R_eff - 1)/β
- e* = i*(γ + μ)/σ
- r* = 1 - s* - e* - i*

where μ is birth/death rate, β is transmission rate, σ is progression rate,
and γ is recovery rate.

# Examples
```julia
# Calculate equilibrium for target disease
target = TargetDiseaseDynamicsParameters(
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2
)
spec = DynamicsParameterSpecification(target)

result = calculate_endemic_equilibrium_proportions(spec, 0.8)

if Try.isok(result)
    props = Try.unwrap(result)
    println("S: \$(props.s_prop), E: \$(props.e_prop), I: \$(props.i_prop), R: \$(props.r_prop)")
else
    println("Error: \$(Try.unwrap_err(result))")
end

# Handle case where R_eff ≤ 1
result = calculate_endemic_equilibrium_proportions(spec, 0.95)  # High vaccination
if Try.iserr(result)
    println("No endemic equilibrium: \$(Try.unwrap_err(result))")
    # Use default initial conditions instead
end
```

# See Also
- [`DynamicsParameters`](@ref): Dynamics parameter type
- [`DynamicsParameterSpecification`](@ref): Specification type
"""
function calculate_endemic_equilibrium_proportions(
        dynamics_params::Union{DynamicsParameters, DynamicsParameterSpecification},
        vaccination_coverage::Float64,
    )
    return calculate_endemic_equilibrium_proportions(
        dynamics_params.R_0,
        dynamics_params.beta_mean,
        dynamics_params.sigma,
        dynamics_params.gamma,
        dynamics_params.mu,
        vaccination_coverage,
    )
end

"""
    calculate_endemic_equilibrium_proportions(R_0, beta_mean, sigma, gamma, mu, vaccination_coverage)

Calculate endemic equilibrium proportions from individual parameters.

This is the low-level implementation that performs the actual calculation.
Most users should use the version that accepts DynamicsParameters instead.

# Arguments
- `R_0::Float64`: Basic reproduction number
- `beta_mean::Float64`: Mean transmission rate
- `sigma::Float64`: Progression rate from E to I
- `gamma::Float64`: Recovery rate
- `mu::Float64`: Birth/death rate
- `vaccination_coverage::Float64`: Vaccination coverage (0-1)

# Returns
- `Try.Ok((s_prop, e_prop, i_prop, r_prop))` on success
- `Try.Err(message)` on failure

# Validation
- Vaccination coverage must be in [0, 1]
- R_eff = R₀(1 - ρ) must be > 1 for endemic equilibrium
- Calculated r_prop must be non-negative

# Examples
```julia
result = calculate_endemic_equilibrium_proportions(
    16.0,  # R_0
    0.0001,  # beta_mean
    0.1,  # sigma (1/10 days)
    0.125,  # gamma (1/8 days)
    4.4e-5,  # mu
    0.8  # vaccination_coverage
)

if Try.isok(result)
    props = Try.unwrap(result)
    # Use equilibrium proportions
else
    # Handle error
    error_msg = Try.unwrap_err(result)
end
```
"""
function calculate_endemic_equilibrium_proportions(
        R_0::Float64,
        beta_mean::Float64,
        sigma::Float64,
        gamma::Float64,
        mu::Float64,
        vaccination_coverage::Float64,
    )
    if vaccination_coverage < 0.0 || vaccination_coverage > 1.0
        return Try.Err(
            "Vaccination coverage must be in [0, 1]. Got $vaccination_coverage"
        )
    end

    R_eff = R_0 * (1.0 - vaccination_coverage)

    if R_eff <= 1.0
        return Try.Err(
            "Effective R₀ must be > 1 for endemic equilibrium. " *
                "Got R_eff = R₀(1 - ρ) = $(round(R_eff; digits = 2))"
        )
    end

    s_prop = 1.0 / R_0
    i_prop = mu * (R_eff - 1.0) / beta_mean
    e_prop = i_prop * (gamma + mu) / sigma
    r_prop = 1.0 - s_prop - e_prop - i_prop

    if r_prop < 0.0
        return Try.Err(
            "Calculated negative recovered proportion. Check parameter consistency."
        )
    end

    return Try.Ok((s_prop = s_prop, e_prop = e_prop, i_prop = i_prop, r_prop = r_prop))
end
