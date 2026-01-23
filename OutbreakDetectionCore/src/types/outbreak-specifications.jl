export OutbreakSpecification

"""
    OutbreakSpecification

Specification for outbreak detection thresholds.

# Fields
- `outbreak_threshold::Integer`: Threshold for outbreak detection
- `minimum_outbreak_duration::Integer`: Minimum duration to be considered outbreak
- `minimum_outbreak_size::Integer`: Minimum size to be considered outbreak
- `dirpath::AbstractString`: Directory path for output

# Constructor
    OutbreakSpecification(outbreak_threshold, minimum_outbreak_duration,
                         minimum_outbreak_size)

Creates an `OutbreakSpecification` with automatically generated directory path.
"""
AutoHashEquals.@auto_hash_equals struct OutbreakSpecification{T1 <: Integer, T2 <: AbstractString}
    outbreak_threshold::T1
    minimum_outbreak_duration::T1
    minimum_outbreak_size::T1
    dirpath::T2
end

function OutbreakSpecification(
        outbreak_threshold,
        minimum_outbreak_duration,
        minimum_outbreak_size
    )
    dirpath = joinpath(
        "min_outbreak_dur_$(minimum_outbreak_duration)",
        "min_outbreak_size_$(minimum_outbreak_size)",
        "outbreak_threshold_$(outbreak_threshold)",
    )

    return OutbreakSpecification(
        outbreak_threshold,
        minimum_outbreak_duration,
        minimum_outbreak_size,
        dirpath,
    )
end
