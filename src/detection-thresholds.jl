# module DetectionThresholds
#
# export create_inc_infec_arr, create_inc_infec_arr!, calculate_outbreak_thresholds

using ProgressMeter
using FLoops
using StatsBase
using UnPack
using LoopVectorization

function create_inc_infec_arr(
    ensemble_inc_vecs, outbreak_specification::OutbreakSpecification
)
    ensemble_inc_arr = zeros(
        Int64, size(ensemble_inc_vecs, 1), 2, size(ensemble_inc_vecs, 2)
    )

    ensemble_thresholds_vec = Vector{Array{Int64,2}}(
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

        abovethresholdrle = rle(@view(ensemble_inc_arr[:, 2, sim]))

        outbreak_thresholds = calculate_outbreak_thresholds(
            abovethresholdrle
        )

        ensemble_thresholds_vec[sim] = classify_all_outbreaks(
            @view(ensemble_inc_arr[:, :, sim]),
            outbreak_thresholds,
            minoutbreakdur,
            minoutbreaksize,
        )
    end
    return nothing
end

function calculate_outbreak_thresholds(outbreakrle)
    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    upperbound_indices = findall(isequal(1), outbreakrle[1])

    outbreak_thresholds = zeros(Int64, length(upperbound_indices), 4)

    @inbounds outbreak_thresholds[:, 1] .= @view(
        outbreakaccum[upperbound_indices]
    )
    map!(
        x -> x - 1 == 0 ? 1 : outbreakaccum[x - 1] + 1,
        outbreak_thresholds[:, 2],
        upperbound_indices,
    )

    return outbreak_thresholds
end

function classify_all_outbreaks(
    incidence_arr,
    all_thresholds_arr,
    minoutbreakdur,
    minoutbreaksize
)
    for (row, (lower, upper)) in pairs(eachrow(all_thresholds_arr[:, 1:2]))
        calculate_period_sum!(
            @view(all_thresholds_arr[row, 3]),
            @view(incidence_arr[lower:upper, 1]),
        )

        all_thresholds_arr[row, 4] = classify_outbreak(
            all_thresholds_arr[row, 3],
            upper,
            lower,
            minoutbreakdur,
            minoutbreaksize,
        )
    end

    return @view(
        all_thresholds_arr[(all_thresholds_arr[:, 1] .!= 0), :]
    )
end

function calculate_period_sum!(outvec, incvec)
    @inbounds outvec .= sum(incvec)
    return nothing
end

function classify_outbreak(
    periodsumvec,
    upper_time,
    lower_time,
    minoutbreakdur,
    minoutbreaksize
)
    if upper_time - lower_time >= minoutbreakdur &&
        periodsumvec >= minoutbreaksize
        return 1
    end
    return 0
end
# end
