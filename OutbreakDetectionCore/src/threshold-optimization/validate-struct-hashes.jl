"""
    validate_struct_hashes_and_get_results(
        data::Dict,
        source_description::String
    ) -> Try.Result{StructVector{OptimizationResult}}

Validate struct hashes in loaded data and return results or prompt user for action.

This function checks if the loaded JLD2 data contains struct version hashes and
validates them against the current struct definitions. It handles three cases:

1. **Hashes match**: Returns the results
2. **Hashes don't match**: Prompts user to confirm re-optimization
3. **No hashes (old format)**: Prompts user to continue, force re-optimization, or quit

# Arguments
- `data::Dict`: Loaded JLD2 data containing "optimization_results" and optionally "struct_hashes"
- `source_description::String`: Description of the data source (e.g., "results file", "checkpoint")
  for user-facing messages

# Returns
- `Try.Ok(results)`: If validation passes and user confirms to continue
- `Try.Ok(StructVector{OptimizationResult}[])`: If user requests re-optimization
- `Try.Err(message)`: If validation fails or user declines to proceed

# Example
```julia
data = JLD2.load(filepath)
if haskey(data, "optimization_results")
    result = validate_struct_hashes_and_get_results(data, "results file")
    if Try.isok(result)
        results = Try.unwrap(result)
        # Use results...
    else
        # Handle error...
    end
end
```

# See Also
- [`get_optimization_struct_hashes`](@ref): Get current struct hashes
- [`confirm_struct_change_reoptimization`](@ref): Prompt user about struct changes
- [`load_previous_optimization_results_structvector`](@ref): Uses this function
- [`load_checkpoint_results_structvector`](@ref): Uses this function
"""
function validate_struct_hashes_and_get_results(
        data::Dict,
        source_description::String
    )
    if !haskey(data, "optimization_results")
        return Try.Err("Data does not contain 'optimization_results' key")
    end

    results = data["optimization_results"]

    # Check if results are the correct type (not ReconstructedMutable from struct changes)
    if !(results isa StructVector{OptimizationResult})
        @warn "Loaded $source_description contains incompatible struct definitions."
        @warn "Results are of type $(typeof(results).name.name) instead of StructVector{OptimizationResult}."
        @warn "This indicates the OptimizationResult struct has changed since results were saved."
        println()
        print("Force re-optimization of all scenarios? (y/N): ")
        response = lowercase(strip(readline()))

        if response in ["y", "yes"]
            @info "User confirmed re-optimization. Invalidating all existing results from $source_description."
            return Try.Ok(StructVector(OptimizationResult[]))
        else
            return Try.Err(
                "User declined re-optimization after incompatible struct type detected in $source_description"
            )
        end
    end

    # Check struct hashes if they exist
    if haskey(data, "struct_hashes")
        stored_hashes = data["struct_hashes"]
        current_hashes = get_optimization_struct_hashes()

        if stored_hashes != current_hashes
            # Struct definitions have changed
            n_existing = length(results)

            if confirm_struct_change_reoptimization(
                    n_existing,
                    stored_hashes,
                    current_hashes
                )
                @info "User confirmed re-optimization. Invalidating all existing results from $source_description."
                return Try.Ok(StructVector(OptimizationResult[]))
            else
                return Try.Err(
                    "User declined re-optimization after struct changes detected in $source_description"
                )
            end
        end
    else
        # Old format without struct hashes - warn user
        @warn "Loaded $source_description does not contain struct version hashes (old format)."
        @warn "Cannot verify struct compatibility. Proceeding with caution."
        println()
        print("Continue using these results, force re-optimization of all scenarios, or quit? (continue/force/quit): ")
        response = lowercase(strip(readline()))

        if response in ["quit", "q", "n", "no"]
            return Try.Err("User declined to use $source_description without struct hashes")
        elseif response in ["force", "f", "reoptimize", "optimization"]
            @info "User requested re-optimization. Invalidating all existing results from $source_description."
            return Try.Ok(StructVector(OptimizationResult[]))
        elseif response in ["continue", "c", "y", "yes"]
            @info "Continuing with existing results from $source_description (struct compatibility not verified)."
            # Fall through to return results below
        else
            @warn "Invalid response '$response'. Defaulting to 'quit' for safety."
            return Try.Err("Invalid user response - defaulting to quit")
        end
    end

    return Try.Ok(results)
end
