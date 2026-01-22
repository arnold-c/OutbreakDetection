export recreate_noise_vecs


"""
    recreate_noise_vecs(noise_spec, ensemble_spec, base_dynamics; verbose=false, seed=1234)

Recreate noise vectors with updated vaccination coverage.

This function generates noise simulation results using the specified dynamical
noise parameters. It uses endemic equilibrium to set initial states for noise
simulations, falling back to simple proportions if equilibrium doesn't exist.

# Arguments
- `noise_spec::DynamicalNoiseSpecification`: Noise specification with vaccination coverage
- `ensemble_spec::EnsembleSpecification`: Ensemble parameters
- `base_dynamics::DynamicsParameterSpecification`: Base dynamics

# Keyword Arguments
- `verbose::Bool`: Print warnings (default: false)
- `seed::Int`: Random seed (default: 1234)

# Returns
- `DynamicalNoiseRun`: Noise simulation results with incidence, Reff, and noise statistics

# Implementation Details
The function:
1. Creates noise-specific dynamics parameters
2. Calculates endemic equilibrium for initial states
3. Falls back to default initial conditions if R_eff ≤ 1
4. Runs ensemble simulation with noise dynamics
5. Calculates mean dynamical noise from results
6. Adds Poisson component if specified
7. Returns DynamicalNoiseRun with all results

# Examples
```julia
# Create noise specification
noise_spec = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 1.0,
    vaccination_coverage = 0.65
)

# Run noise simulation
noise_result = recreate_noise_vecs(
    noise_spec,
    ensemble_spec,
    base_dynamics;
    verbose = true
)

# Access results
println("Mean total noise: \$(noise_result.mean_noise)")
println("Mean dynamical noise: \$(noise_result.mean_dynamic_noise)")
println("Mean Poisson noise: \$(noise_result.mean_poisson_noise)")
```

# See Also
- [`DynamicalNoiseSpecification`](@ref): Noise specification type
- [`calculate_endemic_equilibrium_proportions`](@ref): Endemic equilibrium calculation
- [`create_noise_dynamics_parameters`](@ref): Noise dynamics creation
"""
function recreate_noise_vecs(
        ensemble_specification::EnsembleSpecification,
        vaccination_coverage::Float64;
        verbose = false,
        seed = 1234,
    )
    # Create final EnsembleSpecification with optimal parameters for verification
    UnPack.@unpack state_parameters = ensemble_specification
    UnPack.@unpack init_states, init_state_props = state_parameters
    UnPack.@unpack N = init_states


    UnPack.@unpack state_parameters, time_parameters, nsims = ensemble_spec
    N = state_parameters.init_states.N

    # Create noise-specific dynamics
    noise_dynamics = recreate_noise_dynamics_spec(noise_spec, ensemble_spec)

    # Calculate endemic equilibrium for initial states
    endemic_result = calculate_endemic_equilibrium_proportions(
        noise_dynamics,
        noise_spec.vaccination_coverage
    )

    # Handle endemic equilibrium calculation
    endemic_props = if Try.isok(endemic_result)
        Try.unwrap(endemic_result)
    else
        if verbose
            @warn Try.unwrap_err(endemic_result) *
                "\nDefaulting to no initial infections, s_prop = 1 - vaccination_coverage"
        end
        (
            s_prop = 1.0 - noise_spec.vaccination_coverage,
            e_prop = 0.0,
            i_prop = 0.0,
            r_prop = noise_spec.vaccination_coverage,
        )
    end

    # Create updated state parameters with endemic equilibrium
    updated_state_parameters = StateParameters(;
        N = N,
        s_prop = endemic_props.s_prop,
        e_prop = endemic_props.e_prop,
        i_prop = endemic_props.i_prop,
    )

    updated_ensemble_specification = EnsembleSpecification(
        ensemble_specification.label,
        updated_state_parameters,
        ensemble_specification.time_parameters,
        updated_dynamics_parameter_specification,
        updated_dynamics_parameter_specification,
        ensemble_specification.dynamical_noise_params,
        ensemble_specification.nsims,
        ensemble_specification.dirpath
    )

    updated_dynamics_parameters = DynamicsParameters(;
        beta_mean = updated_dynamics_parameter_specification.beta_mean,
        beta_force = updated_dynamics_parameter_specification.beta_force,
        seasonality = updated_dynamics_parameter_specification.seasonality,
        sigma = updated_dynamics_parameter_specification.sigma,
        gamma = updated_dynamics_parameter_specification.gamma,
        mu = updated_dynamics_parameter_specification.mu,
        annual_births_per_k = updated_dynamics_parameter_specification.annual_births_per_k,
        epsilon = updated_dynamics_parameter_specification.epsilon,
        R_0 = updated_dynamics_parameter_specification.R_0,
        vaccination_coverage = vaccination_coverage
    )

    noise_result = create_noise_vecs(
        updated_dynamical_noise_spec,
        updated_ensemble_specification,
        updated_dynamics_parameters;
        seed = seed
    )

    return noise_result

end
