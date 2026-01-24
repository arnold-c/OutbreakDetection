export save_optimization_results

"""
    save_results_structvector(results, filepath)

Save StructVector results directly to JLD2 file with struct version hashes.

This function saves optimization results along with hash values of the
OptimizationScenario and OptimizationResult struct definitions. These hashes
enable detection of struct changes when loading results, ensuring that cached
results are invalidated if the struct definitions have been modified.
"""
function save_optimization_results(
        results::StructVector{OptimizationResult},
        filepath::String
    )
    # Ensure directory exists
    dir = dirname(filepath)
    !isdir(dir) && mkpath(dir)

    # Create temporary file with .jld2 extension
    temp_filepath = filepath * ".tmp.jld2"

    return try
        # Get current struct hashes
        struct_hashes = get_optimization_struct_hashes()

        # Save StructVector and struct hashes to JLD2
        JLD2.jldsave(
            temp_filepath;
            optimization_results = results,
            struct_hashes = struct_hashes
        )

        # Atomic rename
        mv(temp_filepath, filepath; force = true)
        @info "Saved optimization results to $filepath"

    catch e
        # Clean up temp file if something went wrong
        isfile(temp_filepath) && rm(temp_filepath; force = true)
        rethrow(e)
    end
end
