export optimize_dynamic_noise_params

function optimize_dynamic_noise_params_wrapper(
        ensemble_specification::EnsembleSpecification,
        seir_results::StructVector{SEIRRun},
        target_scaling::Float64,
        optimization_params::NoiseVaccinationOptimizationParameters = NoiseVaccinationOptimizationParameters();
        verbose = false,
        seed = 1234
    )

    # Calculate trimmed mean incidence up to endpoints
    overall_mean = calculate_mean_incidence(seir_results)

    # Call the original optimization function with the computed mean
    return optimize_dynamic_noise_params(
        ensemble_specification,
        overall_mean,
        target_scaling,
        optimization_params;
        verbose = verbose,
        seed = seed
    )
end


"""
    optimize_dynamic_noise_params(target_scaling, mean_target_incidence, noise_spec, ensemble_spec, base_dynamics, optimization_params=NoiseVaccinationOptimizationParameters(); verbose=false, seed=1234)

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
- `optimization_params::NoiseVaccinationOptimizationParameters`: Optimization settings (optional)

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
- Achieved noise differs from target by more than `optimization_params.atol`

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
optimization_params = NoiseVaccinationOptimizationParameters(
    n_sobol_points = 10,
    maxeval = 100,
    atol = 0.5
)
result = optimize_dynamic_noise_params(
    7.0, 3.0, noise_spec, ensemble_spec, base_dynamics, optimization_params;
    verbose = true
)
```

# See Also
- [`NoiseVaccinationOptimizationParameters`](@ref): Optimization configuration
- [`DynamicalNoiseSpecification`](@ref): Noise specification with bounds
- [`calculate_mean_dynamical_noise`](@ref): Objective function
"""
function optimize_dynamic_noise_params(
        ensemble_specification::EnsembleSpecification,
        enddates_vec::Vector{Int64},
        measles_daily_incidence::Float64,
        target_scaling::Float64,
        optimization_params::NoiseVaccinationOptimizationParameters = NoiseVaccinationOptimizationParameters();
        verbose = false,
        seed = 1234
    )
    UnPack.@unpack vaccination_bounds = ensemble_specification.dynamical_noise_specification

    target_noise = target_scaling * mean_target_incidence

    if verbose
        println("Target noise level: $target_noise")
        println("Vaccination bounds: $vaccination_bounds")
        println("Starting multistart optimization with $n_sobol_points points...")
    end

    # Pre-extract ensemble specification components to avoid repeated unpacking
    N = ensemble_specification.state_parameters.init_states.N

    # Define objective function: minimize squared difference from target
    objective = let target_noise = target_noise,
            ensemble_specification = ensemble_specification,
            N = N,
            verbose = verbose
        function (params)
            vaccination_coverage = params[1]

            # Ensure susceptible proportion is valid and will result in positive compartments
            # Use stricter bounds to prevent numerical issues with very small populations
            min_safe_prop = max(1.0 / N, 0.001)  # At least 1 person or 0.1%, whichever is larger
            max_safe_prop = min(1.0 - 1.0 / N, 0.999)  # At most N-1 people or 99.9%, whichever is smaller

            noise_level = calculate_mean_dynamical_noise(
                ensemble_specification,
                vaccination_coverage;
                verbose = verbose,
                seed = seed
            )

            # Return squared error from target
            return (noise_level - target_noise)^2
        end
    end


    # Setup multistart optimization
    problem = MultistartOptimization.MinimizationProblem(
        objective,
        [vaccination_bounds[1]],  # lower bounds
        [vaccination_bounds[2]],   # upper bounds
    )

    # Configure local method
    local_method = MultistartOptimization.NLopt_local_method(
        optimization_params.local_algorithm;
        xtol_rel = optimization_params.xtol_rel,
        xtol_abs = optimization_params.xtol_abs,
        maxeval = optimization_params.maxeval
    )

    # Configure multistart method
    multistart_method = MultistartOptimization.TikTak(optimization_params.n_sobol_points)

    if verbose
        println("Running multistart optimization...")
    end

    # Run optimization
    result = MultistartOptimization.multistart_minimization(
        multistart_method,
        local_method,
        problem
    )

    if !in(result.ret, [:SUCCESS, :XTOL_REACHED, :FTOL_REACHED, :STOPVAL_REACHED]) ||
            sqrt(result.value) > optimization_params.atol

        optimal_vaccination_coverage = result.location[1]

        optimized_noise = calculate_mean_dynamical_noise(
            dynamical_noise_spec,
            optimal_vaccination_coverage,
            ensemble_specification;
            verbose = false,
            seed = seed
        )

        error_msg = "Unsuccessful optimization." *
            "\nReturn code: $(result.ret)" *
            "\nAbsolute difference: $(sqrt(result.value))" *
            "\nTarget value: $(target_noise)" *
            "\nOptimized value: $(optimized_noise)" *
            "\nOptimization parameters: Vax. prop = $(optimal_vaccination_coverage)"


        error(error_msg)
    end

    return result
end
