export recreate_noise_vecs, calculate_mean_dynamical_noise

"""
    recreate_noise_vecs(noise_spec, ensemble_spec, base_dynamics; verbose=false, seed=1234)

Recreate noise vectors with updated vaccination coverage.

This function generates noise simulation results using the specified dynamical
noise parameters. It uses endemic equilibrium to set initial states for noise
simulations, falling back to simple proportions if equilibrium doesn't exist.

# Arguments
- `noise_spec::DynamicalNoise`: Noise specification with vaccination coverage
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
noise_spec = DynamicalNoise(
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
- [`DynamicalNoise`](@ref): Noise specification type
- [`calculate_endemic_equilibrium_proportions`](@ref): Endemic equilibrium calculation
- [`create_noise_dynamics_parameters`](@ref): Noise dynamics creation
"""
function recreate_noise_vecs(
        noise_spec::DynamicalNoise,
        ensemble_spec::EnsembleSpecification,
        base_dynamics::DynamicsParameterSpecification;
        verbose::Bool = false,
        seed::Int = 1234,
    )
    UnPack.@unpack state_parameters, time_parameters, nsims = ensemble_spec
    N = state_parameters.init_states.N

    # Create noise-specific dynamics
    noise_dynamics = create_noise_dynamics_parameters(noise_spec, base_dynamics, N)

    # Calculate endemic equilibrium for initial states
    endemic_result = calculate_endemic_equilibrium_proportions(
        noise_dynamics, noise_spec.vaccination_coverage
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
    noise_state_params = StateParameters(
        N = N,
        s_prop = endemic_props.s_prop,
        e_prop = endemic_props.e_prop,
        i_prop = endemic_props.i_prop,
    )

    # Run ensemble simulation with noise dynamics
    # Note: This assumes run_jump_prob or similar function exists
    # For now, we'll create a placeholder that returns mock data
    # This will need to be updated when integrating with actual simulation code
    noise_seir_results = _run_noise_ensemble_simulation(
        noise_state_params, noise_dynamics, time_parameters, nsims; seed = seed
    )

    # Calculate mean dynamical noise from SEIR results
    mean_dynamic_noise = calculate_mean_incidence(noise_seir_results)

    # Add Poisson component
    mean_poisson_noise = noise_spec.poisson_component * mean_dynamic_noise
    total_mean_noise = mean_dynamic_noise + mean_poisson_noise

    # Extract incidence and Reff vectors
    incidence_vecs = [run.incidence for run in noise_seir_results]
    reff_vecs = [run.Reff for run in noise_seir_results]

    return DynamicalNoiseRun(
        incidence = incidence_vecs,
        Reff = reff_vecs,
        mean_noise = total_mean_noise,
        mean_poisson_noise = mean_poisson_noise,
        mean_dynamic_noise = mean_dynamic_noise,
    )
end

"""
    calculate_mean_dynamical_noise(noise_spec, ensemble_spec, base_dynamics; verbose=false, seed=1234)

Calculate mean dynamical noise for given vaccination coverage.

This is a wrapper function for optimization objectives. It runs noise
simulations and returns only the mean noise level.

# Arguments
- `noise_spec::DynamicalNoise`: Noise specification
- `ensemble_spec::EnsembleSpecification`: Ensemble parameters
- `base_dynamics::DynamicsParameterSpecification`: Base dynamics

# Keyword Arguments
- `verbose::Bool`: Print warnings (default: false)
- `seed::Int`: Random seed (default: 1234)

# Returns
- `Float64`: Mean noise level (dynamical + Poisson components)

# Examples
```julia
# Use in optimization objective
function objective(vaccination_coverage)
    noise = DynamicalNoise(spec, vaccination_coverage)
    noise_level = calculate_mean_dynamical_noise(
        noise,
        ensemble_spec,
        base_dynamics
    )
    return (noise_level - target_noise)^2
end
```

# See Also
- [`recreate_noise_vecs`](@ref): Full noise recreation with all results
- [`optimize_dynamic_noise_params`](@ref): Optimization wrapper
"""
function calculate_mean_dynamical_noise(
        noise_spec::DynamicalNoise,
        ensemble_spec::EnsembleSpecification,
        base_dynamics::DynamicsParameterSpecification;
        verbose::Bool = false,
        seed::Int = 1234,
    )
    result = recreate_noise_vecs(
        noise_spec, ensemble_spec, base_dynamics; verbose = verbose, seed = seed
    )
    return result.mean_noise
end

# Placeholder function for ensemble simulation
# This will be replaced with actual integration to existing simulation code
"""
    _run_noise_ensemble_simulation(state_params, dynamics_params, time_params, nsims; seed=1234)

Placeholder for ensemble simulation integration.

This function will be replaced with proper integration to the existing
ensemble simulation infrastructure. For now, it returns a StructVector
of mock SEIRRun results.

# Note
This is a temporary implementation to allow the noise optimization
infrastructure to be built. It should be replaced with calls to the
actual SEIR simulation code.
"""
function _run_noise_ensemble_simulation(
        state_params::StateParameters,
        dynamics_params::DynamicsParameters,
        time_params::SimTimeParameters,
        nsims::Int;
        seed::Int = 1234,
    )
    # TODO: Replace with actual ensemble simulation
    # For now, return mock data structure
    @warn "Using placeholder ensemble simulation. Replace with actual SEIR simulation."

    # Create mock results with correct structure
    tlength = time_params.tlength
    mock_runs = [
        SEIRRun(
                states = [
                    StaticArrays.SVector{5, Int64}(
                        state_params.init_states.S,
                        state_params.init_states.E,
                        state_params.init_states.I,
                        state_params.init_states.R,
                        state_params.init_states.N,
                    ) for _ in 1:tlength
                ],
                incidence = rand(1:10, tlength),
                Reff = rand(0.5:0.01:2.0, tlength),
            ) for _ in 1:nsims
    ]

    return StructVector{SEIRRun}(mock_runs)
end
