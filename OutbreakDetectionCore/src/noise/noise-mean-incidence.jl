export calculate_mean_incidence

"""
    calculate_mean_incidence(seir_results)

Calculate mean incidence from SEIR simulation results.

This function computes the mean daily incidence across all simulations
and time points. It works efficiently with both StructVector and regular
Vector of SEIRRun results.

# Arguments
- `seir_results::Union{StructVector{SEIRRun}, Vector{SEIRRun}}`: SEIR simulation results

# Returns
- `Float64`: Mean daily incidence across all simulations and time points

# Implementation Notes
When using StructVector, this function takes advantage of efficient column
access. `seir_results.incidence` directly accesses the incidence column as
a `Vector{Vector{Int64}}`, which is more efficient than iterating over
individual SEIRRun structs.

# Examples
```julia
# With StructVector (efficient)
ensemble_results = StructVector{SEIRRun}([run1, run2, run3, ...])
mean_inc = calculate_mean_incidence(ensemble_results)

# With regular Vector (fallback)
ensemble_results = [run1, run2, run3, ...]
mean_inc = calculate_mean_incidence(ensemble_results)

# Use in noise optimization
target_mean = calculate_mean_incidence(target_disease_results)
target_noise = 7.0 * target_mean  # High noise scenario
```

# See Also
- [`SEIRRun`](@ref): Single simulation result type
- [`calculate_mean_dynamical_noise`](@ref): Calculate noise from dynamical simulations
"""
function calculate_mean_incidence(seir_results::StructVector{SEIRRun})
    # StructVector provides efficient column access
    # seir_results.incidence is a Vector{Vector{Int64}}
    all_incidence = vcat(seir_results.incidence...)
    return Statistics.mean(all_incidence)
end

function calculate_mean_incidence(seir_results::Vector{SEIRRun})
    # Fallback for regular Vector
    all_incidence = vcat([result.incidence for result in seir_results]...)
    return Statistics.mean(all_incidence)
end
