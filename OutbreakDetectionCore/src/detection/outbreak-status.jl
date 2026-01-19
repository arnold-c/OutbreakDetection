export create_outbreak_status_vecs

"""
    create_outbreak_status_vecs(
        ensemble_simulation::StructVector{SEIRRun},
        outbreak_specification::OutbreakSpecification
    ) -> (Vector{Vector{Int64}}, Vector{Matrix{Int64}})

Create binary outbreak status vectors and outbreak bounds for each simulation.

For each simulation, determines which time points are during a true outbreak
based on the outbreak specification criteria (threshold, minimum duration,
minimum size).

# Arguments
- `ensemble_simulation`: StructVector of SEIRRun containing incidence vectors
- `outbreak_specification`: OutbreakSpecification with threshold, min duration, min size

# Returns
- `outbreak_status_vecs`: Vector of binary vectors (0/1) indicating outbreak status
- `outbreak_bounds_vecs`: Vector of matrices with outbreak bounds [start, end, duration, size]
"""
function create_outbreak_status_vecs(
        ensemble_simulation::StructVector{SEIRRun},
        outbreak_specification::OutbreakSpecification,
    )
    nsims = length(ensemble_simulation)
    outbreak_status_vecs = Vector{Vector{Int64}}(undef, nsims)
    outbreak_bounds_vecs = Vector{Matrix{Int64}}(undef, nsims)

    for sim in eachindex(ensemble_simulation)
        inc_vec = ensemble_simulation.incidence[sim]

        # Determine which days are above threshold
        abovethreshold_vec = vec(
            inc_vec .>= outbreak_specification.outbreak_threshold
        )

        # Run-length encoding to find consecutive periods
        abovethresholdrle = StatsBase.rle(abovethreshold_vec)

        # Calculate outbreak bounds
        all_outbreak_bounds = calculate_outbreak_thresholds(
            abovethresholdrle; ncols = 5
        )

        # Classify each period as outbreak or not
        outbreak_status = zeros(Int64, length(inc_vec))
        classify_all_outbreaks!(
            outbreak_status,
            all_outbreak_bounds,
            inc_vec,
            outbreak_specification,
        )

        outbreak_status_vecs[sim] = outbreak_status
        outbreak_bounds_vecs[sim] = filter_only_outbreaks(all_outbreak_bounds)
    end

    return outbreak_status_vecs, outbreak_bounds_vecs
end
