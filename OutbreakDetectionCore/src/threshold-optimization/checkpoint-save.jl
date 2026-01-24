"""
    save_optimization_checkpoint(results, checkpoint_dir, batch_idx)

Save checkpoint file atomically for StructVector results with struct version hashes.

This function saves checkpoint results along with hash values of the
OptimizationScenario and OptimizationResult struct definitions to ensure
consistency when loading checkpoints.
"""
function save_optimization_checkpoint(
        results::StructVector{OptimizationResult},
        checkpoint_dir::String,
        checkpoint_output_filename_base::String,
        batch_idx::Int64
    )
    if !isdir(checkpoint_dir)
        mkpath(checkpoint_dir)
    end

    checkpoint_file = joinpath(checkpoint_dir, "$(checkpoint_output_filename_base)$(batch_idx).jld2")
    temp_file = checkpoint_file * ".tmp"

    # Get current struct hashes
    struct_hashes = get_optimization_struct_hashes()

    # save to temp and then move to force an overwrite if back with the same name
    # already exists (though it shouldn't as it uses the datetime in the file name)
    JLD2.jldsave(
        temp_file;
        optimization_results = results,
        batch_idx = batch_idx,
        struct_hashes = struct_hashes
    )
    return mv(temp_file, checkpoint_file; force = true)
end
