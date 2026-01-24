export struct_definition_hash,
    get_optimization_struct_hashes

"""
    struct_definition_hash(::Type{T}) where T -> UInt64

Compute a stable hash of a struct's name, field names, and field types.

This function generates a hash value based on the complete structure definition
(struct name, field names, and their types) of a given struct type. The hash can
be used to detect when a struct definition has changed, which is useful for
invalidating cached results that depend on the struct's structure.

# Arguments
- `::Type{T}`: The struct type to hash

# Returns
- `UInt64`: A hash value representing the struct's definition

# Example
```julia
# Compute hash for OptimizationResult
hash_v1 = struct_definition_hash(OptimizationResult)

# After adding a new field to OptimizationResult...
hash_v2 = struct_definition_hash(OptimizationResult)

# The hashes will differ, indicating a structural change
hash_v1 != hash_v2  # true
```

# Notes
- The hash is based on struct name, field names, and field types, not field values
- Changes to struct name, field order, addition/removal of fields, or type changes
  will all produce different hashes
- Two different structs with identical fields will have different hashes due to
  different struct names
- The hash is stable across Julia sessions for the same struct definition

# See Also
- [`get_optimization_struct_hashes`](@ref): Get hashes for both optimization structs
"""
function struct_definition_hash(::Type{T}) where {T}
    struct_name = nameof(T)
    fields = fieldnames(T)
    types = fieldtypes(T)
    return hash((struct_name, fields, types))
end

"""
    get_optimization_struct_hashes() -> NamedTuple

Get hash values for both OptimizationScenario and OptimizationResult structs.

This function computes and returns the current hash values for the two main
optimization-related structs. These hashes can be stored alongside optimization
results and checked later to detect if the struct definitions have changed.

# Returns
- `NamedTuple`: A named tuple with keys:
  - `scenario_hash::UInt64`: Hash of OptimizationScenario struct definition
  - `result_hash::UInt64`: Hash of OptimizationResult struct definition

# Example
```julia
# Get current struct hashes
current_hashes = get_optimization_struct_hashes()

# Save alongside optimization results
JLD2.jldsave(
    filepath;
    optimization_results = results,
    struct_hashes = current_hashes
)

# Later, check if structs have changed
loaded_hashes = JLD2.load(filepath, "struct_hashes")
if loaded_hashes != get_optimization_struct_hashes()
    @warn "Struct definitions have changed since results were saved"
end
```

# See Also
- [`struct_definition_hash`](@ref): Compute hash for a single struct type
- [`save_optimization_results`](@ref): Saves results with struct hashes
- [`load_previous_optimization_results_structvector`](@ref): Checks struct hashes when loading
"""
function get_optimization_struct_hashes()
    return (
        scenario_hash = struct_definition_hash(OptimizationScenario),
        result_hash = struct_definition_hash(OptimizationResult),
    )
end
