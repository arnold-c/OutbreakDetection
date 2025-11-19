export run_scenario_optimizations, run_missing_scenario_optimizations!,
    run_scenario_optimization, _optimize_single_scenario, _load_result_cache, _save_results

"""
    run_scenario_optimizations(ensemble_specifications, outbreak_specifications, 
                               noise_specifications, outbreak_detection_specifications, 
                               individual_test_specifications, optim_method, 
                               accuracy_functions; kwargs...)

Run optimization for multiple scenarios.

Returns optimization DataFrame or nothing based on parameters.
"""
function run_scenario_optimizations(
        ensemble_specifications,
        outbreak_specifications,
        noise_specifications,
        outbreak_detection_specifications,
        individual_test_specifications,
        optim_method::TMethod = MSO,
        accuracy_functions = [arithmetic_mean, calculate_f_beta_score];
        seed = 1234,
        executor = FLoops.SequentialEx(),
        optimization_filename_base = "alert-threshold-optimization.jld2",
        filedir = outdir("ensemble", "scenario-optimizations"),
        optimization_output_filepath = joinpath(
            filedir,
            string(Dates.now()) * "_" * optimization_filename_base,
        ),
        force = false,
        return_df = true,
        filter_df_results = false,
        save_df = true,
        disable_time_check = false,
        time_per_run_s = 45,
        kwargs...,
    ) where {TMethod <: Type{<:OptimizationMethods}}
    if length(save_df + return_df) == 0
        error(
            "At least one of `save_df` and `return_df` must be `true`. Instead, got `save_df = ` $save_df, and `return_df = ` $return_df"
        )
    end

    if !isdir(filedir)
        mkpath(filedir)
    end

    load_filepath = get_most_recent_optimization_filepath(
        optimization_filename_base,
        filedir,
    )

    optim_df = nothing

    if Try.isok(load_filepath) && !force
        optim_df = JLD2.load(Try.unwrap(load_filepath))["threshold_optim_df"]
    else
        optim_df = DataFrames.DataFrame(
            (
                ensemble_spec = typeof(ensemble_specifications[1])[],
                outbreak_spec = typeof(outbreak_specifications[1])[],
                noise_spec = Union{
                    PoissonNoiseSpecification, DynamicalNoiseSpecification,
                }[],
                outbreak_detection_spec = typeof(
                    outbreak_detection_specifications[1]
                )[],
                test_spec = typeof(individual_test_specifications[1])[],
                optimal_threshold = Float64[],
                optimal_accuracy = Float64[],
                optimization_method = OptimizationMethods[],
                accuracy_function = Function[],
                OT_chars = StructVector{<:OutbreakThresholdChars}[],
            )
        )
    end

    combinations_to_run = calculate_combination_to_run(
        ensemble_specifications,
        outbreak_specifications,
        noise_specifications,
        outbreak_detection_specifications,
        individual_test_specifications,
        optim_method,
        accuracy_functions,
    )

    missing_optimizations = check_missing_scenario_optimizations(
        optim_df,
        combinations_to_run;
        disable_time_check = disable_time_check,
        time_per_run_s = time_per_run_s,
    )

    if Try.iserr(missing_optimizations)
        println(Try.unwrap_err(missing_optimizations))
        if return_df
            if filter_df_results
                println(
                    "Returning previously computed threshold_optim_df, filtered for selected combinations"
                )
                return filter_optim_results(
                    optim_df,
                    combinations_to_run,
                )
            end
            println("Returning previously computed threshold_optim_df")
            return optim_df
        end
        return nothing
    end

    missing_optims_df = Try.unwrap(missing_optimizations)

    if DataFrames.nrow(missing_optims_df) == 0
        @info "🟨 Previous optimal ews_df is the same as the current. No need to save a new version. 🟨"
    end

    run_missing_scenario_optimizations!(
        optim_df,
        missing_optims_df;
        seed = seed,
        executor = executor,
        kwargs...,
    )

    if save_df
        @tagsave(
            optimization_output_filepath,
            Dict("threshold_optim_df" => optim_df)
        )
        @info "🟢 Saved optimal ews_df to $(optimization_output_filepath) 🟢"
    end

    if return_df
        if filter_df_results
            return filter_optim_results(
                optim_df,
                combinations_to_run,
            )
        end
        return optim_df
    end

    return nothing
end

"""
    run_missing_scenario_optimizations!(optim_df, missing_optims_df; 
                                        seed=1234, executor=..., kwargs...)

Run optimizations for missing scenarios and update optim_df in place.
"""
function run_missing_scenario_optimizations!(
        optim_df,
        missing_optims_df;
        seed = 1234,
        executor = FLoops.SequentialEx(),
        kwargs...,
    )
    prog = Progress(DataFrames.nrow(missing_optims_df))

    incidence_sim_grouped_df = DataFrames.groupby(
        missing_optims_df,
        [:ensemble_spec, :outbreak_spec],
    )

    lk = ReentrantLock()
    for inc_gp in incidence_sim_grouped_df
        ensemble_spec = inc_gp[1, :ensemble_spec]
        outbreak_spec = inc_gp[1, :outbreak_spec]

        base_param_dict = @dict(
            ensemble_spec = ensemble_spec,
            outbreak_spec = outbreak_spec,
            seed = seed,
        )

        ensemble_inc_arr, thresholds_vec = setup_optimization(
            base_param_dict
        )

        unique_noise_specs = unique(inc_gp[:, :noise_spec])
        FLoops.@floop executor for noise_spec in unique_noise_specs
            noise_gp = filter(:noise_spec => x -> x == noise_spec, inc_gp)

            println(
                styled"{green:\n=================================================================}"
            )
            println(
                styled"Noise type: {green,inverse: $(getdirpath(noise_spec))}"
            )

            noise_array, noise_means = create_noise_arr(
                noise_spec,
                ensemble_inc_arr;
                ensemble_specification = ensemble_spec,
                seed = seed,
            )

            for detect_test_gp in DataFrames.groupby(
                    noise_gp, [:outbreak_detection_spec, :test_spec]
                )
                outbreak_detection_spec = detect_test_gp[
                    1, :outbreak_detection_spec,
                ]
                individual_test_spec = detect_test_gp[1, :test_spec]

                println(
                    styled"\t\t-> Test specification: {blue: $(get_test_description(individual_test_spec))}, Percent tested: {red,inverse: $(outbreak_detection_spec.percent_tested)}"
                )

                for accuracy_function_df in DataFrames.groupby(
                        detect_test_gp, [:accuracy_function]
                    )
                    accuracy_function = accuracy_function_df[
                        1, :accuracy_function,
                    ]
                    println(
                        styled"\t\t\t\t-> Accuracy Function: {yellow,inverse: $(accuracy_function)}"
                    )

                    obj_inputs = (;
                        ensemble_inc_arr,
                        noise_array,
                        outbreak_detection_specification = outbreak_detection_spec,
                        individual_test_specification = individual_test_spec,
                        thresholds_vec,
                        accuracy_function,
                    )

                    objective_function_closure =
                        x -> objective_function(x, obj_inputs)

                    @assert length(
                        unique(detect_test_gp[:, :optimization_method])
                    ) == 1

                    optim_method = detect_test_gp[1, :optimization_method]

                    optim_minimizer, optim_minimum = optimization_wrapper(
                        objective_function_closure,
                        optim_method;
                        kwargs...,
                    )

                    optimal_outbreak_detection_spec = OutbreakDetectionSpecification(
                        optim_minimizer,
                        outbreak_detection_spec.moving_average_lag,
                        outbreak_detection_spec.percent_visit_clinic,
                        outbreak_detection_spec.percent_clinic_tested,
                        outbreak_detection_spec.alert_method.method_name,
                    )

                    testarr = create_testing_arrs(
                        ensemble_inc_arr,
                        noise_array,
                        optimal_outbreak_detection_spec,
                        individual_test_spec,
                    )[1]

                    OT_chars = calculate_OutbreakThresholdChars(
                        testarr, ensemble_inc_arr, thresholds_vec, noise_means
                    )

                    lock(lk)
                    push!(
                        optim_df,
                        (
                            ensemble_spec,
                            outbreak_spec,
                            noise_spec,
                            outbreak_detection_spec,
                            individual_test_spec,
                            optim_minimizer,
                            1 - optim_minimum,
                            optim_method,
                            accuracy_function,
                            OT_chars,
                        ),
                    )
                    unlock(lk)

                    next!(prog)
                    println("")
                end
            end
        end
    end

    return nothing
end

function filter_optim_results(
        optim_df,
        combinations_to_run;
        scenario_parameter_symbols = [
            :ensemble_spec,
            :outbreak_spec,
            :noise_spec,
            :outbreak_detection_spec,
            :test_spec,
            :optimization_method,
            :accuracy_function,
        ],
    )
    return DataFrames.innerjoin(
        optim_df,
        combinations_to_run;
        on = scenario_parameter_symbols,
    )
end

# ============================================================================
# New optimization wrapper using ScenarioParameters and OptimizationResult
# ============================================================================

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

        # Update progress
        if verbose && !isnothing(progress)
            ProgressMeter.next!(progress)
        end
    end

    # Convert to StructVector for efficient storage
    results_sv = StructVector{OptimizationResult}(results)

    # Save results
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
        # For now, use midpoint of bounds
        # TODO: Implement proper noise level optimization
        vax_coverage = Statistics.mean(scenario.noise_spec.vaccination_bounds)
        DynamicalNoise(scenario.noise_spec, vax_coverage)
    else
        scenario.noise_spec
    end

    # 3. Run ensemble simulation
    # TODO: Integrate with actual ensemble simulation
    # For now, return placeholder

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
