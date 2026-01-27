export confirm_struct_change_reoptimization

"""
    confirm_struct_change_reoptimization(
        n_existing::Int,
        stored_hashes::NamedTuple,
        current_hashes::NamedTuple;
        default_action::Union{Symbol, Nothing} = nothing
    ) -> Symbol

Prompt user for action when struct definitions have changed.

This function is called when loading optimization results reveals that the
OptimizationScenario or OptimizationResult struct definitions have changed
since the results were saved. It displays detailed information about what
changed and asks the user to choose how to proceed.

# Arguments
- `n_existing::Int`: Number of existing results that will be invalidated
- `stored_hashes::NamedTuple`: Hash values stored with the results
- `current_hashes::NamedTuple`: Current hash values of struct definitions

# Keyword Arguments
- `default_action::Union{Symbol, Nothing}`: Default action to take without prompting.
  Must be one of `:continue`, `:force`, or `:quit`. If `nothing` (default), prompts user.

# Returns
- `Symbol`: User's choice (or default_action if provided), one of:
  - `:continue` - Use existing results despite struct changes (risky)
  - `:force` - Force re-optimization of all scenarios
  - `:quit` - Exit without proceeding

# Example
```julia
stored = (scenario_hash = 0x123..., result_hash = 0x456...)
current = get_optimization_struct_hashes()

if stored != current
    # Prompt user
    action = confirm_struct_change_reoptimization(100, stored, current)
    
    # Or use default action
    action = confirm_struct_change_reoptimization(
        100, stored, current;
        default_action = :continue
    )
    
    if action == :force
        # Re-run all optimizations
    elseif action == :continue
        # Use existing results (risky)
    else
        # Exit without running
    end
end
```

# See Also
- [`get_optimization_struct_hashes`](@ref): Get current struct hashes
- [`load_previous_optimization_results_structvector`](@ref): Checks hashes when loading
- [`proceed_with_optimization`](@ref): Similar confirmation for time estimates
"""
function confirm_struct_change_reoptimization(
        n_existing::Int,
        stored_hashes::NamedTuple,
        current_hashes::NamedTuple;
        default_action::Union{Symbol, Nothing} = nothing
    )
    # If default_action is provided, validate and use it
    if !isnothing(default_action)
        if default_action ∉ [:continue, :force, :quit]
            error("default_action must be one of :continue, :force, or :quit, got: $default_action")
        end
        @info "Using default action: $default_action (struct hash mismatch detected)"
        return default_action
    end

    println(StyledStrings.styled"{red:⚠ WARNING: Struct definition changes detected}")
    println()

    # Determine which structs changed
    scenario_changed = stored_hashes.scenario_hash != current_hashes.scenario_hash
    result_changed = stored_hashes.result_hash != current_hashes.result_hash

    if scenario_changed
        println(
            StyledStrings.styled"{yellow:OptimizationScenario} struct has changed:"
        )
        println(
            StyledStrings.styled"  Stored hash:  {cyan:$(stored_hashes.scenario_hash)}"
        )
        println(
            StyledStrings.styled"  Current hash: {cyan:$(current_hashes.scenario_hash)}"
        )
        println()
    end

    if result_changed
        println(StyledStrings.styled"{yellow:OptimizationResult} struct has changed:")
        println(
            StyledStrings.styled"  Stored hash:  {cyan:$(stored_hashes.result_hash)}"
        )
        println(
            StyledStrings.styled"  Current hash: {cyan:$(current_hashes.result_hash)}"
        )
        println()
    end

    println(
        StyledStrings.styled"This means the {cyan:$n_existing} existing results {yellow:may be incompatible}"
    )
    println(
        StyledStrings.styled"with the current struct definitions."
    )
    println()
    println(
        StyledStrings.styled"{bold:Possible causes:} Added/removed fields, changed field types, or reordered fields"
    )
    println()

    print(StyledStrings.styled"{bold:Continue using these results, force re-optimization of all scenarios, or quit?} (continue/force/quit): ")
    response = lowercase(strip(readline()))

    if response in ["quit", "q", "n", "no"]
        return :quit
    elseif response in ["force", "f", "reoptimize", "optimization"]
        return :force
    elseif response in ["continue", "c", "y", "yes"]
        return :continue
    else
        @warn "Invalid response '$response'. Defaulting to 'quit' for safety."
        return :quit
    end
end
