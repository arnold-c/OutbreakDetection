# module DiagTestingFunctions
#
# export create_testing_arr, create_testing_arr!, calculate_tested!,
#     calculate_pos, calculate_pos!, calculate_movingavg, calculate_movingavg!,
#     detectoutbreak, detectoutbreak!, calculate_ot_characterstics,
#     calculate_noutbreaks, calculate_OutbreakThresholdChars,
#     run_OutbreakThresholdChars_creation, OutbreakThresholdChars_creation

using DrWatson
using StatsBase
using FreqTables
using ThreadsX
using FLoops
using StructArrays

# include("detection-thresholds.jl")
# # using .DetectionThresholds
#
# include("structs.jl")
# using .ODStructs

function create_testing_arrs(
    incarr,
    noisearr,
    outbreak_detect_spec::OutbreakDetectionSpecification,
    individual_test_spec::IndividualTestSpecification,
)
    testarr = zeros(Int64, size(incarr, 1), 8, size(incarr, 3))
    testpos_vec = Vector{TestPositivity}(undef, size(incarr, 3))
    ntested_worker_vec = Vector{Int64}(undef, size(incarr, 1))

    create_testing_arrs!(
        testarr,
        testpos_vec,
        ntested_worker_vec,
        incarr,
        noisearr,
        outbreak_detect_spec.detection_threshold,
        outbreak_detect_spec.moving_average_lag,
        outbreak_detect_spec.percent_tested,
        outbreak_detect_spec.test_result_lag,
        individual_test_spec.sensitivity,
        individual_test_spec.specificity,
    )

    return testarr, StructArray(testpos_vec)
end

function create_testing_arrs!(
    testarr,
    testpos_vec,
    ntested_worker_vec,
    incarr,
    noisearr,
    detectthreshold,
    moveavglag,
    perc_tested,
    testlag,
    testsens,
    testspec,
)
    ntested = size(testarr, 1)

    for sim in axes(incarr, 3)
        # Number of infectious individuals tested
        calculate_tested!(
            @view(testarr[:, 1, sim]), @view(incarr[:, 1, sim]), perc_tested
        )

        # Number of noise individuals tested
        calculate_tested!(
            @view(testarr[:, 2, sim]), @view(noisearr[:, 1, sim]), perc_tested
        )

        # Number of TOTAL individuals tested
        @turbo @. @views ntested_worker_vec .=
            testarr[:, 1, sim] + testarr[:, 2, sim]

        # Number of test positive INFECTED individuals
        calculate_true_positives!(
            @view(testarr[:, 3, sim]),
            @view(testarr[:, 1, sim]),
            ntested,
            testlag,
            testsens,
        )

        # Number of test positive NOISE individuals
        calculate_noise_positives!(
            @view(testarr[:, 4, sim]),
            @view(testarr[:, 2, sim]),
            ntested,
            testlag,
            testspec,
        )

        # Number of test positive TOTAL individuals
        @. testarr[:, 5, sim] =
            @view(testarr[:, 3, sim]) + @view(testarr[:, 4, sim])

        # Calculate moving average of TOTAL test positives
        calculate_movingavg!(
            @view(testarr[:, 6, sim]),
            @view(testarr[:, 5, sim]),
            testlag, moveavglag;
            Float = false,
        )

        # TOTAL Test positive individuals trigger outbreak response
        detectoutbreak!(
            @view(testarr[:, 7, sim]),
            @view(testarr[:, 5, sim]),
            @view(testarr[:, 6, sim]),
            detectthreshold,
        )

        # Triggered outbreak equal to actual outbreak status
        @. testarr[:, 8, sim] =
            @view(testarr[:, 7, sim]) == @view(incarr[:, 4, sim])

        # Posterior prob of infectious / total test tests performed
        testpos_vec[sim] = TestPositivity(
            @view(testarr[:, 5, sim]),
            ntested_worker_vec,
            @view(testarr[:, 7, sim])
        )
    end

    return nothing
end

function calculate_tested!(outvec, invec, perc_tested)
    @. outvec = round(invec * perc_tested)
end

function calculate_positives(tested_vec, lag, tested_multiplier)
    ntested = length(tested_vec)
    npos = zeros(Int64, ntested)

    calculate_positives!(
        npos,
        tested_vec,
        ntested,
        lag,
        tested_multiplier,
    )

    return npos
end

function calculate_noise_positives!(outvec, tested_vec, tlength, lag, spec)
    tested_multiplier = 1.0 - spec
    calculate_positives!(outvec, tested_vec, tlength, lag, tested_multiplier)
    return nothing
end

function calculate_true_positives!(outvec, tested_vec, tlength, lag, sens)
    calculate_positives!(outvec, tested_vec, tlength, lag, sens)
    return nothing
end

function calculate_positives!(
    npos_vec, tested_vec, ntested, lag, tested_multiplier
)
    for day in eachindex(tested_vec)
        if day + lag <= ntested
            npos_vec[day + lag] = Int64(
                round(tested_vec[day] * tested_multiplier)
            )
        end
    end
    return nothing
end

function calculate_movingavg(invec, testlag, avglag)
    outvec = zeros(Float64, size(invec, 1), 1)

    calculate_movingavg!(outvec, invec, testlag, avglag)

    return outvec
end

function calculate_movingavg!(outvec, invec, testlag, avglag; Float = true)
    if Float
        avgfunc =
            (invec, day, avglag) -> mean(@view(invec[(day - avglag + 1):day]))
    else
        avgfunc =
            (invec, day, avglag) -> Int64(round(
                mean(@view(invec[(day - avglag + 1):day]))
            ))
    end
    for day in eachindex(invec)
        if day >= testlag + avglag + 1
            outvec[day] = avgfunc(invec, day, avglag)
        end
    end
end

function detectoutbreak(incvec, avgvec, threshold)
    outbreak = zeros(Int64, length(incvec))

    detectoutbreak!(outbreak, incvec, avgvec, threshold)

    return outbreak
end

function detectoutbreak!(outbreakvec, incvec, avgvec, threshold)
    @. outbreakvec = ifelse(incvec >= threshold || avgvec >= threshold, 1, 0)

    return nothing
end

function calculate_test_positivity(
    true_positive_vec, total_positive_vec, detection_vec, agg_days
)
    @views outvec = zeros(Float64, length(true_positive_vec) ÷ agg_days, 2)
    @inbounds for i in axes(outvec, 1)
        start_ind = 1 + (i - 1) * agg_days
        end_ind = start_ind + (agg_days - 1)

        @views total_positive_sum = sum(total_positive_vec[start_ind:end_ind])
        @views true_positive_sum = sum(true_positive_vec[start_ind:end_ind])
        @views num_outbreak_days = sum(detection_vec[start_ind:end_ind])
        agg_outbreak_status = num_outbreak_days >= agg_days / 2 ? 1 : 0

        outvec[i, 1] = true_positive_sum / total_positive_sum
        outvec[i, 2] = agg_outbreak_status
    end
    return outvec
end

function calculate_ot_characterstics(testvec, infecvec)
    crosstab = freqtable(testvec, infecvec)

    tp = freqtable_error_default_zero(crosstab, 1, 1)
    fp = freqtable_error_default_zero(crosstab, 1, 0)
    tn = freqtable_error_default_zero(crosstab, 0, 0)
    fn = freqtable_error_default_zero(crosstab, 0, 1)

    sens = tp / (tp + fn)
    spec = tn / (tn + fp)

    ppv = tp / (tp + fp)
    npv = tn / (tn + fn)

    return crosstab, tp, tn, fp, fn, sens, spec, ppv, npv
end

# TODO: write tests that check the function works when crosstab produces a 2x2 table,
# a 2x1, and a 1x1 table
function freqtable_error_default_zero() end

function freqtable_error_default_zero(
    freqtable, testing_class::T, actual_class::T
) where {T<:Integer}
    return try
        freqtable[Name(testing_class), Name(actual_class)]
    catch
        0
    end
end

function calculate_noutbreaks(outbreakrle)
    return length(findall(==(1), outbreakrle[1]))
end

function calculate_OutbreakThresholdChars(testarr, infecarr)
    OT_chars = ThreadsX.map(axes(infecarr, 3)) do sim
        outbreakrle = rle(@view(infecarr[:, 4, sim]))
        detectrle = rle(@view(testarr[:, 7, sim]))

        OutbreakThresholdChars(
            calculate_ot_characterstics(
                @view(testarr[:, 7, sim]), @view(infecarr[:, 4, sim])
            )...,
            calculate_noutbreaks(outbreakrle),
            calculate_noutbreaks(detectrle),
            reduce(hcat, collect(calculate_outbreak_thresholds(outbreakrle))),
            reduce(hcat, collect(calculate_outbreak_thresholds(detectrle))),
        )
    end

    return StructArray(OT_chars)
end

# end
