using StatsBase: StatsBase

export create_inc_infec_arr, create_inc_infec_arr!,
    calculate_outbreak_thresholds

"""
    create_inc_infec_arr(ensemble_inc_vecs, outbreak_specification)

Create incidence and infection arrays with outbreak classification.

Returns ensemble_inc_arr and ensemble_thresholds_vec.
"""
function create_inc_infec_arr(
        ensemble_inc_vecs, outbreak_specification::OutbreakSpecification
    )
    ensemble_inc_arr = zeros(
        Int64, size(ensemble_inc_vecs, 1), 3, size(ensemble_inc_vecs, 2)
    )

    ensemble_thresholds_vec = Vector{Array{Int64, 2}}(
        undef, size(ensemble_inc_vecs, 2)
    )

    create_inc_infec_arr!(
        ensemble_inc_arr,
        ensemble_thresholds_vec,
        ensemble_inc_vecs,
        outbreak_specification.outbreak_threshold,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )

    return ensemble_inc_arr, ensemble_thresholds_vec
end

"""
    create_inc_infec_arr!(ensemble_inc_arr, ensemble_thresholds_vec, 
                          ensemble_inc_vecs, outbreakthreshold, 
                          minoutbreakdur, minoutbreaksize)

In-place version of create_inc_infec_arr.
"""
function create_inc_infec_arr!(
        ensemble_inc_arr, ensemble_thresholds_vec, ensemble_inc_vecs,
        outbreakthreshold, minoutbreakdur,
        minoutbreaksize,
    )
    @inbounds for sim in axes(ensemble_inc_vecs, 2)
        convert_svec_to_matrix!(
            @view(ensemble_inc_arr[:, 1, sim]),
            @view(ensemble_inc_vecs[:, sim])
        )

        ensemble_inc_arr[:, 2, sim] .=
            @view(ensemble_inc_arr[:, 1, sim]) .>= outbreakthreshold

        abovethresholdrle = StatsBase.rle(@view(ensemble_inc_arr[:, 2, sim]))

        outbreak_thresholds = calculate_outbreak_thresholds(
            abovethresholdrle; ncols = 5
        )

        classify_all_outbreaks!(
            @view(ensemble_inc_arr[:, 3, sim]),
            outbreak_thresholds,
            @view(ensemble_inc_arr[:, 1, sim]),
            minoutbreakdur,
            minoutbreaksize,
        )

        ensemble_thresholds_vec[sim] = filter_only_outbreaks(
            outbreak_thresholds
        )
    end
    return nothing
end

"""
    calculate_outbreak_thresholds(outbreakrle; ncols=5)

Calculate outbreak threshold bounds from run-length encoding.
"""
function calculate_outbreak_thresholds(outbreakrle; ncols = 5)
    @assert ncols >= 3

    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    upperbound_indices = findall(isequal(1), outbreakrle[1])

    outbreak_thresholds = zeros(Int64, length(upperbound_indices), ncols)

    @inbounds outbreak_thresholds[:, 2] .= @view(
        outbreakaccum[upperbound_indices]
    )
    map!(
        x -> x - 1 == 0 ? 1 : outbreakaccum[x - 1] + 1,
        @view(outbreak_thresholds[:, 1]),
        upperbound_indices,
    )

    return outbreak_thresholds
end
