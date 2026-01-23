export validate_all_simulations_have_outbreaks

"""
    validate_all_simulations_have_outbreaks(
        outbreak_thresholds::StructVector{OutbreakThresholds},
        ensemble_specification::EnsembleSpecification,
        outbreak_specification::OutbreakSpecification
    )

Validate that all simulations in the ensemble have at least one outbreak.

This function checks each simulation in the ensemble to ensure it contains at least
one outbreak period. If any simulation lacks outbreaks, it returns an error wrapped
in a Try.Err containing detailed information about the failure, including the
ensemble specification, outbreak specification, and the specific simulation number
that failed.

# Arguments
- `outbreak_thresholds::StructVector{OutbreakThresholds}`: Outbreak thresholds for each simulation
- `ensemble_specification::EnsembleSpecification`: Ensemble configuration used to generate simulations
- `outbreak_specification::OutbreakSpecification`: Outbreak detection criteria

# Returns
- `Try.Ok{Nothing}`: If all simulations have at least one outbreak
- `Try.Err`: If any simulation lacks outbreaks, containing error message with:
  - Simulation number that failed
  - Ensemble specification details
  - Outbreak specification details

# Example
```julia
outbreak_thresholds = calculate_outbreak_thresholds(ensemble_simulation, outbreak_spec)
result = validate_all_simulations_have_outbreaks(
    outbreak_thresholds,
    ensemble_spec,
    outbreak_spec
)

# Handle the result
if Try.isok(result)
    # Continue with optimization
else
    error(Try.unwrap_err(result))
end
```

# See Also
- [`calculate_outbreak_thresholds`](@ref): Function that generates outbreak thresholds
- [`OutbreakThresholds`](@ref): Type containing outbreak information per simulation
- [`EnsembleSpecification`](@ref): Ensemble configuration type
- [`OutbreakSpecification`](@ref): Outbreak detection criteria type
"""
function validate_all_simulations_have_outbreaks(
        outbreak_thresholds::StructVector{OutbreakThresholds},
        ensemble_specification::EnsembleSpecification,
        outbreak_specification::OutbreakSpecification
    )
    for (sim_idx, outbreak_threshold) in pairs(outbreak_thresholds)
        if isempty(outbreak_threshold.lower_bounds)
            error_msg = """
            Simulation $sim_idx has no outbreaks detected.

            This indicates that the outbreak specification parameters are too strict
            for the ensemble simulation parameters, resulting in at least one simulation
            with no qualifying outbreak periods.

            Ensemble Specification:
            $(ensemble_specification)

            Outbreak Specification:
            $(outbreak_specification)

            Simulation Index: $sim_idx
            Total Simulations: $(length(outbreak_thresholds))

            Please adjust either:
            1. Outbreak detection criteria (lower thresholds, shorter duration, smaller size)
            2. Ensemble parameters (higher transmission, longer simulation time)
            """
            return Try.Err(error_msg)
        end
    end

    return Try.Ok(nothing)
end
