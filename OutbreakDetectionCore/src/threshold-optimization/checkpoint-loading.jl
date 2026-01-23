export load_checkpoint_results_structvector

"""
    load_checkpoint_results_structvector(filedir)

Load the most recent valid checkpoint file from the directory as StructVector.
Falls back to older checkpoints if the most recent one is corrupted.
"""
function load_checkpoint_results_structvector(filedir::String)
    checkpoint_dir = joinpath(filedir, "checkpoints")

    if !isdir(checkpoint_dir)
        @warn "$checkpoint_dir isn't a directory. Returning an empty StructVector to force the recreation of all scenarios."
        print("Continue? (y/N): ")
        response = readline()
        if lowercase(strip(response)) in ["y", "yes"]
            return Try.Ok(StructVector(OptimizationResult[]))
        else
            return Try.Err("Quitting optimization")
        end
    end

    checkpoint_files = filter(
        f -> endswith(f, ".jld2") && contains(f, "checkpoint_"),
        readdir(checkpoint_dir)
    )

    if isempty(checkpoint_files)
        @warn "Failed to load the previous checkpoints in $checkpoint_dir. Returning an empty StructVector to force the recreation of all scenarios."
        print("Continue? (y/N): ")
        response = readline()
        if lowercase(strip(response)) in ["y", "yes"]
            return Try.Ok(StructVector(OptimizationResult[]))
        else
            return Try.Err("Quitting optimization")
        end
    end

    # Sort checkpoint files by modification time (most recent first)
    sorted_files = sort(
        checkpoint_files,
        by = f -> mtime(joinpath(checkpoint_dir, f)),
        rev = true
    )

    for file in sorted_files
        filepath = joinpath(checkpoint_dir, file)
        try
            data = JLD2.load(filepath)
            if haskey(data, "optimization_results")
                results = data["optimization_results"]
                if results isa StructVector{OptimizationResult}
                    if file != sorted_files[1]
                        @info "Loaded checkpoint from $file (most recent checkpoint failed)"
                    end
                    return Try.Ok(results)
                end
            end
        catch e
            @warn "Failed to load checkpoint file $file: $e. Trying next checkpoint..."
        end
    end

    @warn "All checkpoint files failed to load. Returning an empty StructVector to force the recreation of all scenarios."
    print("Continue? (y/N): ")
    response = readline()
    if lowercase(strip(response)) in ["y", "yes"]
        return Try.Ok(StructVector(OptimizationResult[]))
    else
        return Try.Err("Quitting optimization")
    end
    return Try.Err("All checkpoint files failed to load")
end
