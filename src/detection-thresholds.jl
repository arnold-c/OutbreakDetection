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
        Int64, size(ensemble_inc_vecs, 1), 4, size(ensemble_inc_vecs, 2)
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

function create_inc_infec_arr!(
    ensemble_inc_arr, ensemble_thresholds_vec, ensemble_inc_vecs,
    outbreakthreshold, minoutbreakdur,
    minoutbreaksize,
)
    @inbounds for sim in axes(ensemble_inc_vecs, 2)
        all_thresholds_arr = zeros(Int64, size(ensemble_inc_vecs, 1), 4)

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

        @inbounds for (lower, upper) in zip(
            outbreak_thresholds.lowers,
            outbreak_thresholds.uppers
        )
            all_thresholds_arr[lower, 1] = lower
            all_thresholds_arr[lower, 2] = upper

            period_sum = sum(@view(ensemble_inc_arr[lower:upper, 1, sim]))

            all_thresholds_arr[lower, 3] = period_sum

            if upper - lower >= minoutbreakdur &&
                period_sum >= minoutbreaksize
                outbreak_class = 1
            else
                outbreak_class = 0
            end

            all_thresholds_arr[lower, 4] = outbreak_class

            # calculate_period_sum!(
            #     @view(ensemble_inc_arr[lower:upper, 3, sim]),
            #     @view(ensemble_inc_arr[lower:upper, 1, sim])
            # )
            # classify_outbreak!(
            #     @view(ensemble_inc_arr[lower:upper, 4, sim]),
            #     @view(all_thresholds_arr[lower, :]),
            #     ensemble_inc_arr[lower, 3, sim],
            #     upper,
            #     lower,
            #     minoutbreakdur,
            #     minoutbreaksize,
            # )
        end

        ensemble_thresholds_vec[sim] = @view(all_thresholds_arr[
            (all_thresholds_arr[:, 1] .!= 0), :,
        ])

    end
    return nothing
end

function calculate_outbreak_thresholds(outbreakrle)
    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    upperbound_indices = findall(isequal(1), outbreakrle[1])
    @inbounds outbreakuppers = outbreakaccum[upperbound_indices]
    @inbounds outbreaklowers = map(
        x -> x - 1 == 0 ? 1 : outbreakaccum[x - 1] + 1, upperbound_indices
    )

    return (lowers = outbreaklowers, uppers = outbreakuppers)
end

function calculate_period_sum!(outvec, incvec)
    @inbounds outvec .= sum(incvec)
    return nothing
end

function classify_outbreak!(
    outvec, thresholds_vec, periodsumvec, upper_time, lower_time,
    minoutbreakdur,
    minoutbreaksize,
)
    if upper_time - lower_time >= minoutbreakdur &&
        periodsumvec >= minoutbreaksize
        @inbounds outvec .= 1
        thresholds_vec[1] = lower_time
        thresholds_vec[2] = upper_time

        # BUG: Can't use push! as causes a data race issue
        # push!(outbreak_thresholds_vec, (lower = lower_time, upper = upper_time))
    end
    return nothing
end
# end
