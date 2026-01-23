export load_previous_optimization_results_structvector

"""
    load_previous_optimization_results_structvector(
        filedir::String,
        filename_base::String
    ) -> StructVector{OptimizationResult}

Load previously completed optimization results from the most recent file or checkpoint.

This function attempts to load optimization results from a directory, first trying to
load the most recent completed results file, and falling back to checkpoint files if
the completed results cannot be loaded. It handles various failure modes gracefully
by returning an empty `StructVector` when no valid results can be found.

# Arguments
- `filedir::String`: Directory path containing optimization result files
- `filename_base::String`: Base name of the optimization file (without datetime prefix or extension)

# Returns
- `StructVector{OptimizationResult}`: Previously completed optimization results.
  Returns an empty `StructVector` if no valid results are found.

# Algorithm
1. Validates that the directory exists (errors if not)
2. Attempts to find and load the most recent completed results file
3. If no completed results file exists, attempts to load from checkpoint files
4. If loading fails, attempts to load from checkpoint files as fallback
5. Returns loaded results or an empty `StructVector` if all attempts fail

# Example
```julia
# Load previous optimization results
results = load_previous_optimization_results_structvector(
    "/path/to/results",
    "threshold-optimization"
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
- Warns when falling back to checkpoint files

# See Also
- [`get_most_recent_optimization_filepath`](@ref): Finds the most recent optimization file
- [`load_checkpoint_results_structvector`](@ref): Loads results from checkpoint files
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

        # Try to load as StructVector first (new format)
        if haskey(data, "optimization_results")
            return Try.Ok(data["optimization_results"])
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
        return load_checkpoint_results_structvector(filedir)
    end
end
