using DrWatson
@quickactivate "OutbreakDetection"

using StatsBase
using FreqTables
using ThreadsX

includet(srcdir("Julia/DrWatson-helpers.jl"))
include(funsdir("structs.jl"))

function create_testing_arr(
    incarr, noisearr, perc_tested, testlag, testsens, testspec, detectthreshold,
    moveavglag,
)
    testarr = zeros(Int64, size(incarr, 1), 6, size(incarr, 3))
    posoddsarr = zeros(Float64, size(incarr, 1), 2, size(incarr, 3))

    create_testing_arr!(
        testarr, incarr, noisearr, posoddsarr, perc_tested, testlag, testsens,
        testspec,
        detectthreshold, moveavglag,
    )

    return testarr
end

function create_testing_arr!(
    testarr, incarr, noisearr, posoddsarr, perc_tested, testlag, testsens,
    testspec,
    detectthreshold, moveavglag,
)
    ntested = size(testarr, 1)

    # prog = Progress(size(incarr, 3))
    @floop for sim in 1:size(incarr, 3)
        # Number of infectious individuals tested
        calculate_tested!(testarr, 1, incarr, perc_tested, sim)

        # Number of noise individuals tested
        calculate_tested!(testarr, 2, noisearr, perc_tested, sim)

        # Number of test positive INFECTED individuals
        calculate_pos!(
            @view(testarr[:, 3, sim]),
            @view(testarr[:, 1, sim]),
            ntested,
            testlag,
            testsens,
            testspec;
            noise = false,
        )

        # Number of test positive NOISE individuals
        calculate_pos!(
            @view(testarr[:, 4, sim]),
            @view(testarr[:, 2, sim]),
            ntested,
            testlag,
            testsens,
            testspec;
            noise = true,
        )

        # Number of test positive TOTAL individuals
        @. testarr[:, 5, sim] =
            @view(testarr[:, 3, sim]) + @view(testarr[:, 4, sim])

        # Calculate moving average of TOTAL test positives
        calculate_movingavg!(
            @view(testarr[:, 5, sim]),
            @view(testarr[:, 6, sim]),
            testlag, moveavglag;
            Float = false,
        )

        # TOTAL Test positive individuals trigger outbreak response
        detectoutbreak!(
            @view(testarr[:, 7, sim]),
            @view(testarr[:, 5, sim]),
            @view(testarr[:, 6, sim]),
            detectthreshold, moveavglag,
        )

        # # Posterior prob of infectious / total test positive
        @. @view(posoddsarr[:, 1, sim]) =
            @view(testarr[:, 3, sim]) / @view(testarr[:, 5, sim])
        calculate_movingavg!(
            @view(posoddsarr[:, 1, sim]),
            @view(posoddsarr[:, 2, sim]),
            testlag, moveavglag,
        )

        # Triggered outbreak equal to actual outbreak status
        @. testarr[:, 8, sim] =
            @view(testarr[:, 7, sim]) == @view(incarr[:, 4, sim])

        # next!(prog)
    end

    return nothing
end

function calculate_tested!(outarr, outarr_ind, inarr, perc_tested, sim)
    @. outarr[:, outarr_ind, sim] = round(@view(inarr[:, 1, sim]) * perc_tested)
end

function calculate_pos(
    tested_vec,
    lag,
    sens,
    spec;
    noise = false,
)
    ntested = length(tested_vec)
    npos = zeros(Int64, ntested)

    calculate_pos!(
        npos,
        tested_vec,
        ntested,
        lag,
        sens,
        spec;
        noise = noise,
    )

    return npos
end

function calculate_pos!(
    npos_vec,
    tested_vec,
    ntested,
    lag,
    sens,
    spec;
    noise = false
)
    if noise
        for day in eachindex(tested_vec)
            if day + lag <= ntested
                npos_vec[day + lag] = Int64(
                    round(tested_vec[day] * (1.0 - spec))
                )
            end
        end
    else
        for day in eachindex(tested_vec)
            if day + lag <= ntested
                npos_vec[day + lag] = Int64(
                    round(tested_vec[day] * sens)
                )
            end
        end
    end

    return nothing
end

function calculate_movingavg(invec, testlag, avglag)
    outvec = zeros(Float64, size(invec, 1), 1)

    calculate_movingavg!(invec, outvec, testlag, avglag)

    return outvec
end

function calculate_movingavg!(invec, outvec, testlag, avglag; Float = true)
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

function detectoutbreak(incvec, avgvec, threshold, avglag)
    outbreak = zeros(Int64, length(incvec))

    detectoutbreak!(outbreak, incvec, avgvec, threshold, avglag)

    return outbreak
end

function detectoutbreak!(outbreakvec, incvec, avgvec, threshold)
    @. outbreakvec = ifelse(incvec >= threshold || avgvec >= threshold, 1, 0)

    return nothing
end

function calculate_ot_characterstics(testarr, infecarr, ind)
    crosstab = freqtable(testarr[:, 5, ind], infecarr[:, 4, ind])

    tp = crosstab[2, 2]
    tn = crosstab[1, 1]
    fp = crosstab[2, 1]
    fn = crosstab[1, 2]

    sens = tp / (tp + fn)
    spec = tn / (tn + fp)

    ppv = tp / (tp + fp)
    npv = tn / (tn + fn)

    return crosstab, tp, tn, fp, fn, sens, spec, ppv, npv
end

function calculate_noutbreaks(outbreakrle)
    return length(findall(==(1), outbreakrle[1]))
end

function calculate_OutbreakThresholdChars(testarr, infecarr)
    OT_chars = ThreadsX.map(axes(infecarr, 3)) do sim
        outbreakrle = rle(@view(infecarr[:, 4, sim]))
        detectrle = rle(@view(testarr[:, 7, sim]))

        OutbreakThresholdChars(
            calculate_ot_characterstics(testarr, infecarr, sim)...,
            calculate_noutbreaks(outbreakrle),
            calculate_noutbreaks(detectrle),
            reduce(hcat, collect(calculate_outbreak_thresholds(outbreakrle))),
            reduce(hcat, collect(calculate_outbreak_thresholds(detectrle))),
        )
    end

    return OT_chars
end
