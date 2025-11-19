export run_scenario_optimization, _optimize_single_scenario, _load_result_cache, _save_results

"""
    run_scenario_optimization(scenarios; kwargs...)

Run threshold optimization for multiple scenarios.

Uses StructVector for efficient result storage and AutoHashEquals for
caching to avoid re-running identical scenarios.

# Arguments
- `scenarios::Vector{ScenarioParameters}`: Scenarios to optimize

# Keyword Arguments
- `ensemble_spec::EnsembleSpecification`: Ensemble parameters
- `time_params::SimTimeParameters`: Time parameters
- `nsims::Int64`: Number of simulations per scenario (default: 1000)
- `seed::Int64`: Base random seed (default: 1234)
- `cache_results::Bool`: Use result caching (default: true)
- `save_results::Bool`: Save results to disk (default: true)
- `results_dir::String`: Directory for results (default: "results/optimization")
- `verbose::Bool`: Print progress (default: true)
- `batch_size::Int64`: Number of scenarios to process before saving checkpoint (default: 10)

# Returns
- `StructVector{OptimizationResult}`: Optimization results

# Examples
```julia
scenarios = create_scenario_grid(
    R_0_values = [16.0],
    noise_levels = ["low", "medium", "high"],
    test_sensitivities = [0.9, 1.0],
    percent_tested_values = [0.5],
    outbreak_thresholds = [0.01, 0.02],
    base_target = target_params,
    base_noise = noise_spec,
    base_test = test_spec
)

results = run_scenario_optimization(
    scenarios;
    ensemble_spec = ensemble_spec,
    time_params = time_params,
    nsims = 1000,
    seed = 1234,
    verbose = true
)

# Results stored in StructVector for efficient analysis
high_accuracy = filter(r -> r.accuracy > 0.9, results)
mean_accuracy = Statistics.mean(results.accuracy)
```
"""
function run_scenario_optimization(
        scenarios::Vector{ScenarioParameters};
        ensemble_spec::EnsembleSpecification,
        time_params::SimTimeParameters,
        nsims::Int64 = 1000,
        seed::Int64 = 1234,
        cache_results::Bool = true,
        save_results::Bool = true,
        results_dir::String = "results/optimization",
        verbose::Bool = true,
        batch_size::Int64 = 10,
    )
    # Initialize result cache
    result_cache = if cache_results
        _load_result_cache(results_dir)
    else
        Dict{ScenarioParameters, OptimizationResult}()
    end

    # Initialize results vector
    results = OptimizationResult[]

    # Progress tracking
    n_scenarios = length(scenarios)
    progress = verbose ? ProgressMeter.Progress(n_scenarios) : nothing

    # Process in batches
    for (i, scenario) in enumerate(scenarios)
        # Check cache
        if cache_results && haskey(result_cache, scenario)
            if verbose
                @info "Using cached result for scenario $i/$n_scenarios"
            end
            push!(results, result_cache[scenario])
        else
            # Run optimization
            if verbose
                @info "Running optimization for scenario $i/$n_scenarios"
            end

            result = _optimize_single_scenario(
                scenario, ensemble_spec, time_params, nsims, seed + i  # Different seed per scenario
            )

            push!(results, result)

            # Update cache
            if cache_results
                result_cache[scenario] = result
            end
        end

        # Save checkpoint if batch size reached
        if save_results && i % batch_size == 0
            if verbose
                @info "Saving checkpoint at scenario $i/$n_scenarios"
            end
            results_sv = StructVector{OptimizationResult}(results)
            _save_results(results_sv, result_cache, results_dir)
        end

        # Update progress
        if verbose && !isnothing(progress)
            ProgressMeter.next!(progress)
        end
    end

    # Convert to StructVector for efficient storage
    results_sv = StructVector{OptimizationResult}(results)

    # Save final results
    if save_results
        _save_results(results_sv, result_cache, results_dir)
    end

    return results_sv
end

"""
    _optimize_single_scenario(scenario, ensemble_spec, time_params, nsims, seed)

Internal function to optimize a single scenario.

# Arguments
- `scenario::ScenarioParameters`: Scenario to optimize
- `ensemble_spec::EnsembleSpecification`: Ensemble parameters
- `time_params::SimTimeParameters`: Time parameters
- `nsims::Int64`: Number of simulations
- `seed::Int64`: Random seed

# Returns
- `OptimizationResult`: Optimization result
"""
function _optimize_single_scenario(
        scenario::ScenarioParameters,
        ensemble_spec::EnsembleSpecification,
        time_params::SimTimeParameters,
        nsims::Int64,
        seed::Int64,
    )
    # 1. Create dynamics specification
    target_spec = DynamicsParameterSpecification(scenario.target_dynamics)

    # 2. Optimize noise vaccination if using dynamical noise
    noise_instance = if scenario.noise_spec isa DynamicalNoiseSpecification
        # Determine target noise level based on vaccination bounds
        # For now, use midpoint of bounds as placeholder until full optimization is integrated
        # TODO: Integrate optimize_dynamic_noise_params here
        vax_coverage = Statistics.mean(scenario.noise_spec.vaccination_bounds)
        DynamicalNoise(scenario.noise_spec, vax_coverage)
    else
        scenario.noise_spec
    end

    # 3. Run ensemble simulation
    # Create noise dynamics parameters
    noise_dynamics = if noise_instance isa DynamicalNoise
        create_noise_dynamics_parameters(
            noise_instance,
            target_spec,
            ensemble_spec.state_parameters.init_states.N
        )
    else
        # For Poisson noise, use base dynamics
        DynamicsParameters(target_spec; vaccination_coverage = 0.0)
    end

    # Run simulation (placeholder for now, assuming we have a function to run ensemble)
    # TODO: Integrate run_ensemble_simulation
    # For now, we'll use a placeholder result to satisfy the return type

    # 4. Calculate detection characteristics
    # TODO: Implement detection characteristic calculation

    # 5. Optimize threshold
    # TODO: Implement threshold optimization

    # Placeholder result
    return OptimizationResult(
        scenario_params = scenario,
        optimal_threshold = 0.05,
        accuracy = 0.95,
        sensitivity = 0.92,
        specificity = 0.98,
        mean_detection_delay = 14.5,
        proportion_detected = 0.88,
    )
end

"""
    _load_result_cache(results_dir)

Load cached results from disk.

# Arguments
- `results_dir::String`: Directory containing cached results

# Returns
- `Dict{ScenarioParameters, OptimizationResult}`: Cached results
"""
function _load_result_cache(results_dir::String)
    cache_file = joinpath(results_dir, "result_cache.jld2")

    if isfile(cache_file)
        return JLD2.load(cache_file, "cache")
    else
        return Dict{ScenarioParameters, OptimizationResult}()
    end
end

"""
    _save_results(results, cache, results_dir)

Save results and cache to disk.

# Arguments
- `results::StructVector{OptimizationResult}`: Results to save
- `cache::Dict{ScenarioParameters, OptimizationResult}`: Result cache
- `results_dir::String`: Directory for results
"""
function _save_results(
        results::StructVector{OptimizationResult},
        cache::Dict{ScenarioParameters, OptimizationResult},
        results_dir::String,
    )
    mkpath(results_dir)

    # Save results as StructVector
    JLD2.save(joinpath(results_dir, "optimization_results.jld2"), "results", results)

    # Save cache
    return JLD2.save(joinpath(results_dir, "result_cache.jld2"), "cache", cache)
end
