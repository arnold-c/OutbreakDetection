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
    @unpack outbreak_threshold, minimum_outbreak_duration, minimum_outbreak_size = outbreak_specification

    incarr = create_inc_infec_arr(
        ensemble_jump_arr, outbreak_threshold, minimum_outbreak_duration, minimum_outbreak_size
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
    incarr, ensemblejumparr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)
    prog = Progress(size(ensemblejumparr, 3))
    @floop for sim in axes(ensemblejumparr, 3)
        # Copy new infections to array
        incarr[:, 1, sim] = @view(ensemblejumparr[1, :, sim])
        # Calculate if new infection is above or below threshold
        incarr[:, 2, sim] = @view(incarr[:, 1, sim]) .>= outbreakthreshold

        # Calculate the total number of infections above threshold in a consecutive string of days
        # Calculate the number of consecutive days of infection above or below threshold
        above5rle = rle(@view(incarr[:, 2, sim]))

        ## Calculate upper and lower indices of consecutive days of infection
        above5lowers, above5uppers = calculate_outbreak_thresholds(above5rle)

        for (lower, upper) in zip(above5lowers, above5uppers)
            # Calculate number of infections between lower and upper indices
            period_sum = sum(@view(incarr[lower:upper, 1, sim]))
            incarr[lower:upper, 3, sim] .= period_sum

            # Determine if there is an outbreak between lower and upper indices
            if upper - lower >= minoutbreakdur && period_sum >= minoutbreaksize
                incarr[lower:upper, 4, sim] .= 1
            end
        end

        next!(prog)
    end
end

function calculate_outbreak_thresholds(outbreakrle)
    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    outbreakuppers = outbreakaccum[findall(==(1), outbreakrle[1])]
    outbreaklowers = filter(
        x -> x <= maximum(outbreakuppers),
        outbreakaccum[findall(==(0), outbreakrle[1])] .+ 1,
    )

    return (outbreaklowers, outbreakuppers)
end

# end