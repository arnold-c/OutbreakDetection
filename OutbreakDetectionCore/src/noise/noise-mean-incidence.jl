export calculate_mean_dynamical_noise

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
        ensemble_specification::EnsembleSpecification,
        mean_vaccination_coverage::Float64;
        verbose = false,
        seed = 1234
    )

    noise_result = recreate_noise_vecs(
        ensemble_specification,
        mean_vaccination_coverage;
        verbose = verbose,
        seed = seed
    )

    return noise_result.mean_noise
end
