# module DetectionThresholds
#
# export create_inc_infec_arr, create_inc_infec_arr!, calculate_outbreak_thresholds

using ProgressMeter
using FLoops
using StatsBase
using UnPack
using LoopVectorization

function create_inc_infec_arr(
    ensemble_jump_arr, outbreak_specification::OutbreakSpecification
)
    ensemble_inc_arr = zeros(
        Int64, size(ensemble_jump_arr, 1), 4, size(ensemble_jump_arr, 3)
    );

    return create_inc_infec_arr!(
        ensemble_inc_arr,
        ensemble_jump_arr,
        outbreak_specification.outbreak_threshold,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )

end

function create_inc_infec_arr!(
    incarr, ensemblejumparr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)
    incarr[:, 1, :] .= @view(ensemblejumparr[:, 1, :])
    for sim in axes(ensemblejumparr, 3)
        incarr[:, 2, sim] .= @view(incarr[:, 1, sim]) .>= outbreakthreshold

        abovethresholdrle = rle(@view(incarr[:, 2, sim]))
        abovethresholdlowers, abovethresholduppers = calculate_outbreak_thresholds(abovethresholdrle)

        for (lower, upper) in zip(abovethresholdlowers, abovethresholduppers)
            calculate_period_sum!(incarr, lower, upper, sim)
            classify_outbreak!(incarr, lower, upper, sim, minoutbreakdur, minoutbreaksize)
        end
    end
    return nothing
end

function calculate_outbreak_thresholds(outbreakrle)
    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    upperbound_indices = findall(isequal(1), outbreakrle[1])
    outbreakuppers = outbreakaccum[upperbound_indices]
    outbreaklowers = map(
        x -> x - 1 == 0 ? 1 : outbreakaccum[x - 1] + 1, upperbound_indices
    )

    return (outbreaklowers, outbreakuppers)
end

function calculate_period_sum!(incarr, lower, upper, sim)
    incarr[lower:upper, 3, sim] .= sum(
        view(incarr, lower:upper, 1, sim)
    )
    return nothing
end

function classify_outbreak!(
    incarr, lower, upper, sim, minoutbreakdur, minoutbreaksize
)
    if upper - lower >= minoutbreakdur && incarr[lower, 3, sim] >= minoutbreaksize

        incarr[lower:upper, 4, sim] .= 1
    end
    return nothing
end
# end
