export optimize_dynamic_noise_params

"""
    optimize_dynamic_noise_params(target_scaling, mean_target_incidence, noise_spec, ensemble_spec, base_dynamics, opt_params=NoiseVaccinationOptimizationParameters(); verbose=false, seed=1234)

Optimize vaccination coverage to achieve target noise level.

This is the primary function for determining vaccination coverage in dynamical
noise simulations. It finds the vaccination coverage that produces a mean noise
level equal to `target_scaling * mean_target_incidence`.

# Workflow
1. User specifies target_scaling (e.g., 7.0 for "high noise")
2. User provides mean_target_incidence from target disease simulations
3. Optimization finds vaccination coverage where:
   mean_noise_incidence = target_scaling × mean_target_incidence

**No burnin or Reff targets needed** - this is purely a noise level matching problem.

# Algorithm
Uses multistart optimization with:
- Sobol sequences for global search initialization
- NLopt local optimization (default: BOBYQA)
- Squared error objective: (achieved_noise - target_noise)²

# Arguments
- `target_scaling::Float64`: Multiplicative factor for target noise (e.g., 1.0 for "low", 7.0 for "high")
- `mean_target_incidence::Float64`: Mean daily incidence from target disease simulations
- `noise_spec::DynamicalNoiseSpecification`: Noise specification with vaccination_bounds
- `ensemble_spec::EnsembleSpecification`: Ensemble parameters for noise simulations
- `base_dynamics::DynamicsParameterSpecification`: Base disease dynamics (not used for Reff)
- `opt_params::NoiseVaccinationOptimizationParameters`: Optimization settings (optional)

# Keyword Arguments
- `verbose::Bool`: Print optimization progress (default: false)
- `seed::Int`: Random seed for reproducibility (default: 1234)

# Returns
Named tuple with:
- `optimal_vaccination::Float64`: Optimal vaccination coverage that achieves target noise
- `mean_noise::Float64`: Achieved mean noise level
- `target_noise::Float64`: Target noise level (target_scaling × mean_target_incidence)
- `difference::Float64`: Absolute difference from target
- `optimization_result`: Full optimization result object

# Errors
Throws an error if:
- Optimization fails to converge
- Achieved noise differs from target by more than `opt_params.atol`

# Examples
```julia
# Define noise specification with vaccination bounds
noise_spec = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 1.0,
    vaccination_bounds = [0.0, 0.8]  # Search range for optimization
)

# Optimize to find vaccination coverage for "high noise" (7x target incidence)
result = optimize_dynamic_noise_params(
    7.0,   # target_scaling: want noise 7x higher than target disease
    3.0,   # mean_target_incidence: from target disease simulations
    noise_spec,
    ensemble_spec,
    base_dynamics
)

# Result: vaccination coverage that produces mean noise ≈ 21.0 (7.0 × 3.0)
println("Use vaccination coverage: \$(result.optimal_vaccination)")
println("This produces mean noise: \$(result.mean_noise)")

# Fast optimization for testing
opt_params = NoiseVaccinationOptimizationParameters(
    n_sobol_points = 10,
    maxeval = 100,
    atol = 0.5
)
result = optimize_dynamic_noise_params(
    7.0, 3.0, noise_spec, ensemble_spec, base_dynamics, opt_params;
    verbose = true
)
```

# See Also
- [`NoiseVaccinationOptimizationParameters`](@ref): Optimization configuration
- [`DynamicalNoiseSpecification`](@ref): Noise specification with bounds
- [`calculate_mean_dynamical_noise`](@ref): Objective function
"""
function optimize_dynamic_noise_params(
        target_scaling::Float64,
        mean_target_incidence::Float64,
        noise_spec::DynamicalNoiseSpecification,
        ensemble_spec::EnsembleSpecification,
        base_dynamics::DynamicsParameterSpecification,
        opt_params::NoiseVaccinationOptimizationParameters = NoiseVaccinationOptimizationParameters();
        verbose::Bool = false,
        seed::Int = 1234,
    )
    UnPack.@unpack vaccination_bounds = noise_spec
    UnPack.@unpack n_sobol_points, local_algorithm, xtol_rel, xtol_abs, maxeval, atol =
        opt_params

    target_noise = target_scaling * mean_target_incidence

    if verbose
        println("Target noise level: $target_noise")
        println("Vaccination bounds: $vaccination_bounds")
        println("Starting multistart optimization with $n_sobol_points points...")
    end

    # Define objective function
    objective = function (params)
        vaccination_coverage = params[1]

        # Create noise instance with this vaccination coverage
        noise_instance = DynamicalNoise(noise_spec, vaccination_coverage)

        # Calculate mean noise
        noise_level = calculate_mean_dynamical_noise(
            noise_instance, ensemble_spec, base_dynamics; seed = seed
        )

        if verbose
            println(
                "  Vax: $(round(vaccination_coverage; digits = 4)), " *
                    "Noise: $(round(noise_level; digits = 2)), " *
                    "Target: $(round(target_noise; digits = 2))",
            )
        end

        # Return squared error
        return (noise_level - target_noise)^2
    end

    # Setup multistart optimization
    problem = MultistartOptimization.MinimizationProblem(
        objective,
        [vaccination_bounds[1]],  # lower bounds
        [vaccination_bounds[2]],   # upper bounds
    )

    # Configure local method
    local_method = MultistartOptimization.NLopt_local_method(
        local_algorithm; xtol_rel = xtol_rel, xtol_abs = xtol_abs, maxeval = maxeval
    )

    # Configure multistart method
    multistart_method = MultistartOptimization.TikTak(n_sobol_points)

    if verbose
        println("Running multistart optimization...")
    end

    # Run optimization
    result = MultistartOptimization.multistart_minimization(
        multistart_method, local_method, problem
    )

    # Extract results
    optimal_vaccination = result.location[1]
    achieved_noise = calculate_mean_dynamical_noise(
        DynamicalNoise(noise_spec, optimal_vaccination),
        ensemble_spec,
        base_dynamics;
        seed = seed,
    )
    difference = abs(achieved_noise - target_noise)

    # Check convergence
    if !in(result.ret, [:SUCCESS, :XTOL_REACHED, :FTOL_REACHED, :STOPVAL_REACHED]) ||
            difference > atol
        error(
            "Optimization failed or did not converge.\n" *
                "Return code: $(result.ret)\n" *
                "Optimal vaccination: $(optimal_vaccination)\n" *
                "Target noise: $(target_noise)\n" *
                "Achieved noise: $(achieved_noise)\n" *
                "Difference: $(difference)\n" *
                "Tolerance: $(atol)",
        )
    end

    if verbose
        println("\nOptimization successful!")
        println("Optimal vaccination coverage: $(round(optimal_vaccination; digits = 4))")
        println("Target noise: $(round(target_noise; digits = 2))")
        println("Achieved noise: $(round(achieved_noise; digits = 2))")
        println("Difference: $(round(difference; digits = 4))")
    end

    return (
        optimal_vaccination = optimal_vaccination,
        mean_noise = achieved_noise,
        target_noise = target_noise,
        difference = difference,
        optimization_result = result,
    )
end
