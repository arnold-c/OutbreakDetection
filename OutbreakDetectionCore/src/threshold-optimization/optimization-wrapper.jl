export run_scenario_threshold_optimization

"""
    gridsearch_structvector(specification_vecs, data_arrs; kwargs...)

Grid search optimization using StructVector for efficient parallel processing.
Reuses functions from multistart optimization where possible.

# Arguments
- `specification_vecs`: Named tuple containing specification vectors including grid parameters
- `data_arrs`: Named tuple containing ensemble data arrays

# Keyword Arguments
- `filedir`: Directory for saving results
- `gridsearch_filename_base`: Base filename for results
- `executor`: Executor for loops (default: `FLoops.ThreadedEx()`)
- `batch_size`: Batch size for processing scenarios (default: 50)
- `force`: Force recomputation even if results exist (default: false)
- `return_results` Return the results (default: true)
- `save_results`: Save results to file (default: true)
- `verbose`: Print progress information (default: true)
- `disable_time_check`: Skip time estimation confirmation (default: false)
"""
function run_scenario_threshold_optimization(
        specification_vecs::ScenarioSpecificationVecs;
        # Optimization control
        dynamic_noise_optimization_parameters::NoiseVaccinationOptimizationParameters = NoiseVaccinationOptimizationParameters(),
        threshold_optimization_parameters::ThresholdOptimizationParameters = ThresholdOptimizationParameters(),
        # File management
        filedir = outdir("ensemble", "threshold-gridsearch"),
        gridsearch_filename_base = "threshold-gridsearch-structvector.jld2",
        gridsearch_output_filepath = joinpath(
            filedir,
            string(Dates.now()) * "_" * gridsearch_filename_base,
        ),
        checkpoint_dir = outdir("ensemble", "threshold-gridsearch", "checkpoints"),
        checkpoint_output_filename_base = joinpath(
            filedir,
            string(Dates.now()) * "_" * "checkpoint_batch_",
        ),
        # Control options
        scheduler = :dynamic, #:serial, :greedy, :static, :dynamic
        force = false,
        return_results = true,
        save_results = true,
        save_checkpoints = true,
        save_checkpoint_num = 1,
        verbose = true,
        disable_time_check = false,
        seconds_per_scenario = 0.025
    )
    if !isdir(filedir)
        mkpath(filedir)
    end
    @assert scheduler in [:dynamic, :static, :greedy, :serial]


    # Create all grid search scenarios as StructVector
    # This includes all combinations with grid parameters
    all_scenarios = create_scenarios_structvector(specification_vecs)

    if verbose
        n_total_scenarios = length(all_scenarios)
        println(StyledStrings.styled"{green:Starting Grid Search with StructVector}")
        println(StyledStrings.styled"Total grid points: {yellow:$n_total_scenarios}")
    end

    # Load existing results (including checkpoints)
    existing_results = if force
        StructVector(OptimizationResult[])
    else
        load_previous_optimization_results_structvector(filedir, gridsearch_filename_base)
    end

    n_existing = length(existing_results)
    if verbose
        println(StyledStrings.styled"Found {cyan:$n_existing} existing results")
    end

    # Find missing scenarios - reuse function from multistart
    missing_scenarios = find_missing_scenarios(all_scenarios, existing_results)
    n_missing = length(missing_scenarios)

    if verbose
        println(StyledStrings.styled"Missing grid points to evaluate: {yellow:$n_missing}")
    end

    # Check with user if needed
    if n_missing > 0
        if !proceed_with_optimization(
                missing_scenarios;
                disable_time_check = disable_time_check,
                seconds_per_scenario = seconds_per_scenario
            )
            @info "Threshold optimization cancelled by user"
            return return_results ? existing_results : nothing
        end
    else
        @info "All optimization scenarios already evaluated"
        return return_results ? existing_results : nothing
    end

    # Run grid search for missing scenarios in batches
    start_time = Dates.now()
    new_results = evaluate_missing_optimizations(
        missing_scenarios;
        save_results = save_results,
        save_checkpoints = save_checkpoints,
        save_checkpoint_num = save_checkpoint_num,
        checkpoint_dir = checkpoint_dir,
        checkpoint_output_filename_base = checkpoint_output_filename_base,
        scheduler = scheduler,
        verbose = verbose,
        dynamic_noise_optimization_parameters = dynamic_noise_optimization_parameters,
        threshold_optimization_parameters = threshold_optimization_parameters,
    )
    end_time = Dates.now()

    # Wrap in try block just in case the datetime arithmetic fails so it doesn't prevent
    # results being saved and returned
    try
        time_taken = Dates.seconds(end_time - start_time)
        println(StyledStrings.styled"Took {yellow:$time_taken} seconds to complete {blue:$n_missing} missing scenarios: {green:$(time_taken/n_missing)} seconds per scenario")
    catch e
    end


    # Combine existing and new results
    BangBang.append!!(existing_results, new_results)

    # Save final results - reuse function from multistart
    if save_results
        save_optimization_results(existing_results, gridsearch_output_filepath)

        # Clean up checkpoints after successful save
        cleanup_checkpoints(checkpoint_dir)
    end

    if verbose
        n_final = length(existing_results)
        println(StyledStrings.styled"{green:Grid search complete! Total results: {yellow:$n_final}}")
    end

    return return_results ? existing_results : nothing
end
