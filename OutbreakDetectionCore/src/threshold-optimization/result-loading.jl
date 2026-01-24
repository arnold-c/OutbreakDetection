export load_previous_optimization_results_structvector

"""
    load_previous_optimization_results_structvector(
        filedir::String,
        filename_base::String
    ) -> Try.Result{StructVector{OptimizationResult}}

Load previously completed optimization results from the most recent file or checkpoint.

This function attempts to load optimization results from a directory, first trying to
load the most recent completed results file, and falling back to checkpoint files if
the completed results cannot be loaded. It also validates that the struct definitions
(OptimizationScenario and OptimizationResult) haven't changed since the results were
saved. If struct changes are detected, it prompts the user to confirm re-optimization.

# Arguments
- `filedir::String`: Directory path containing optimization result files
- `filename_base::String`: Base name of the optimization file (without datetime prefix or extension)
- `checkpoint_dir::String`: Directory path containing checkpoint files

# Returns
- `Try.Result{StructVector{OptimizationResult}}`: Previously completed optimization results.
  Returns an empty `StructVector` if no valid results are found or if struct changes
  are detected and user confirms re-optimization.

# Algorithm
1. Validates that the directory exists (errors if not)
2. Attempts to find and load the most recent completed results file
3. Checks if struct hashes match current struct definitions
4. If struct changes detected, prompts user for confirmation
5. If no completed results file exists or loading fails, attempts to load from checkpoint files
6. Returns loaded results or an empty `StructVector` if all attempts fail

# Example
```julia
# Load previous optimization results
results = Try.unwrap(
    load_previous_optimization_results_structvector(
        "/path/to/results",
        "threshold-optimization",
        "/path/to/checkpoints"
    )
)

# Check if any results were loaded
if isempty(results)
    @info "No previous results found, starting fresh"
else
    @info "Loaded \$(length(results)) previous results"
end
```

# Error Cases
- Throws an error if `filedir` is not a valid directory
- Returns empty `StructVector` if no optimization files are found
- Returns empty `StructVector` if all files fail to load
- Returns empty `StructVector` if struct changes detected and user confirms re-optimization
- Returns `Try.Err` if user declines to proceed with re-optimization
- Warns when falling back to checkpoint files

# See Also
- [`get_most_recent_optimization_filepath`](@ref): Finds the most recent optimization file
- [`load_checkpoint_results_structvector`](@ref): Loads results from checkpoint files
- [`confirm_struct_change_reoptimization`](@ref): Prompts user about struct changes
- [`get_optimization_struct_hashes`](@ref): Gets current struct hashes
- [`OptimizationResult`](@ref): Type of results stored in the returned `StructVector`
"""
function load_previous_optimization_results_structvector(
        filedir::String,
        filename_base::String,
        checkpoint_dir::String
    )
    if !isdir(filedir)
        error("The provided file directory $filedir isn't a valid directory. Check `filedir` again.")
    end

    # Get most recent file - reuse function from multistart
    load_filepath = get_most_recent_optimization_filepath(filename_base, filedir)

    checkpoint_warning_message = "Failed to load previous completed results so attempting to load the most recent checkpoint file"

    if Try.iserr(load_filepath)
        # Try to load checkpoint files
        @warn checkpoint_warning_message
        println(Try.unwrap_err(load_filepath))
        return load_checkpoint_results_structvector(checkpoint_dir)
    end

    try
        data = JLD2.load(Try.unwrap(load_filepath))

        # Validate struct hashes and get results
        if haskey(data, "optimization_results")
            return validate_struct_hashes_and_get_results(data, "results file")
        end

        @warn "Failed to load the previous results in $load_filepath. Returning an empty StructVector to force the recreation of all scenarios."
        print("Continue? (y/N): ")
        response = readline()
        if lowercase(strip(response)) in ["y", "yes"]
            return Try.Ok(StructVector(OptimizationResult[]))
        else
            return Try.Err("Quitting optimization")
        end
    catch
        # Try to load checkpoint files as fallback
        @warn checkpoint_warning_message
        return load_checkpoint_results_structvector(checkpoint_dir)
    end
end
