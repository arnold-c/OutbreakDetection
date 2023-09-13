# module DetectionThresholds
#
# export create_inc_infec_arr, create_inc_infec_arr!, calculate_outbreak_thresholds

using ProgressMeter
using FLoops
using StatsBase
using UnPack

function create_inc_infec_arr(
    ensemble_jump_arr, outbreak_specification::OutbreakSpecification
)
    incarr = create_inc_infec_arr(
        ensemble_jump_arr,
        outbreak_specification.outbreak_threshold,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )

    return incarr
end

function create_inc_infec_arr(
    ensemble_jump_arr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)
    incarr = zeros(
        Int64, size(ensemble_jump_arr, 2), 4, size(ensemble_jump_arr, 3)
    )

    create_inc_infec_arr!(
        incarr,
        ensemble_jump_arr,
        outbreakthreshold,
        minoutbreakdur,
        minoutbreaksize,
    )

    return incarr
end

function create_inc_infec_arr!(
    incarr, ensemble_jump_arr, outbreak_specification::OutbreakSpecification
)
    create_inc_infec_arr!(
        incarr,
        ensemble_jump_arr,
        outbreak_specification.outbreak_threshold,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )
    return nothing
end

function create_inc_infec_arr!(
    incarr, ensemblejumparr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)
    incarr[:, 2, :] .= @view(ensemblejumparr[1, :, :]) .>= outbreakthreshold

    for sim in axes(ensemblejumparr, 3)
        above5rle = rle(@view(incarr[:, 2, sim]))

        ## Calculate upper and lower indices of consecutive days of infection
        above5lowers, above5uppers = calculate_outbreak_thresholds(above5rle)

        for (lower, upper) in zip(above5lowers, above5uppers)
            calculate_period_sum!(incarr, ensemblejumparr, lower, upper, sim)
            classify_outbreak!(
                incarr, lower, upper, sim, minoutbreakdur, minoutbreaksize
            )
        end
    end
end

function calculate_outbreak_thresholds(outbreakrle)
    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    upperbound_indices = findall(isequal(1), outbreakrle[1])
    outbreakuppers = @view(outbreakaccum[upperbound_indices])
    outbreaklowers = map(
        x -> x - 1 == 0 ? 1 : outbreakaccum[x - 1] + 1, upperbound_indices
    )

    return (outbreaklowers, outbreakuppers)
end

function calculate_period_sum!(incarr, jumparr, lower, upper, sim)
    return incarr[lower:upper, 3, sim] .= sum(
        @view(jumparr[1, lower:upper, sim])
    )
end

function classify_outbreak!(
    incarr, lower, upper, sim, minoutbreakdur, minoutbreaksize
)
    if lower - upper >= minoutbreakdur &&
        incarr[lower, 3, sim] >= minoutbreaksize
        incarr[lower:upper, 4, sim] .= 1
    end
end
# end
