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

function create_testing_arr(
    incarr,
    noisearr,
    perc_tested,
    testlag,
    testsens,
    testspec,
    detectthreshold,
    moveavglag,
)
    testarr = zeros(Int64, size(incarr, 1), 8, size(incarr, 3))
    posoddsarr = zeros(Float64, size(incarr, 1), 2, size(incarr, 3))

    create_testing_arr!(
        testarr, incarr, noisearr, posoddsarr, perc_tested, testlag, testsens,
        testspec,
        detectthreshold, moveavglag,
    )

    return testarr
end

function create_testing_arr!(
    testarr,
    incarr,
    noisearr,
    posoddsarr,
    perc_tested,
    testlag,
    testsens,
    testspec,
    detectthreshold,
    moveavglag,
)
    ntested = size(testarr, 1)

    # prog = Progress(size(incarr, 3))
    @floop for sim in axes(incarr, 3)
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
            detectthreshold,
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

function calculate_pos(tested_vec, lag, sens, spec; noise = false)
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
    npos_vec, tested_vec, ntested, lag, sens, spec; noise = false
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

    return StructArray(OT_chars)
end

function run_OutbreakThresholdChars_creation(
    dict_of_OTchars_params; progress = true
)
    if progress
        prog = Progress(length(dict_of_OTchars_params))
    end

    for OTChars_params in dict_of_OTchars_params
        @produce_or_load(
            OutbreakThresholdChars_creation,
            OTChars_params,
            datadir(
                "seasonal-infectivity-import",
                "tau-leaping",
                "N_$(OTChars_params[:N])",
                "r_$(OTChars_params[:init_states_prop][:r_prop])",
                "nsims_$(OTChars_params[:nsims])",
                "births_per_k_$(OTChars_params[:births_per_k])",
                "beta_force_$(OTChars_params[:beta_force])",
                "tmax_$(OTChars_params[:time_p].tmax)",
                "tstep_$(OTChars_params[:time_p].tstep)",
                "noise_$(OTChars_params[:noise_spec].noise_type)",
                "min_outbreak_dur_$(OTChars_params[:outbreak_spec].minimum_outbreak_duration)",
                "min_outbreak_size_$(OTChars_params[:outbreak_spec].minimum_outbreak_size)",
                "outbreak_threshold_$(OTChars_params[:outbreak_spec].outbreak_threshold)",
                "detectthreshold_$(OTChars_params[:outbreak_detect_spec].detection_threshold)",
                "testlag_$(OTChars_params[:outbreak_detect_spec].test_result_lag)",
                "moveavglag_$(OTChars_params[:outbreak_detect_spec].moving_average_lag)",
                "perc_tested_$(OTChars_params[:outbreak_detect_spec].percent_tested)",
                "testsens_$(OTChars_params[:ind_test_spec].sensitivity)",
                "testspec_$(OTChars_params[:ind_test_spec].specificity)",
            );
            prefix = "SEIR_tau_sol",
            filename = savename(
                OTChars_params;
                allowedtypes = (Symbol, Dict, String, Real),
                accesses = [
                    :N,
                    :init_states_prop,
                    :nsims,
                    :time_p,
                    :births_per_k,
                    :beta_force,
                    :noise_spec,
                    :outbreak_spec,
                    :outbreak_detect_spec,
                    :ind_test_spec,
                ],
                expand = [
                    "init_states_prop",
                    "outbreak_spec",
                    "outbreak_detect_spec",
                    "ind_test_spec",
                ],
                sort = false,
            ),
            loadfile = false
        )
        if progress
            next!(prog)
        end
    end
end

function OutbreakThresholdChars_creation(OT_chars_param_dict)
    @unpack N,
    init_states_prop,
    nsims,
    time_p,
    births_per_k,
    beta_force,
    noise_spec,
    outbreak_spec,
    outbreak_detect_spec,
    ind_test_spec = OT_chars_param_dict

    ensemble_spec = EnsembleSpecification(
        ("seasonal-infectivity-import", "tau-leaping"),
        N,
        init_states_prop[:r_prop],
        nsims,
        births_per_k,
        beta_force,
        time_p,
    )

    ensemble_sol = get_ensemble_file(
        "sol", ensemble_spec
    )

    @unpack ensemble_jump_arr = ensemble_sol

    @info "Creating Incidence Array"
    incarr = create_inc_infec_arr(ensemble_jump_arr, outbreak_spec)

    testarr = zeros(Int64, size(incarr, 1), 8, size(incarr, 3))

    posoddsarr = zeros(Float64, size(incarr, 1), 2, size(incarr, 3))

    @unpack noise_array = noise_spec

    @unpack detection_threshold,
    moving_average_lag,
    percent_tested,
    test_result_lag = outbreak_detect_spec

    @unpack sensitivity, specificity = ind_test_spec

    @info "Creating Testing Array"
    create_testing_arr!(
        testarr,
        incarr,
        noise_array,
        posoddsarr,
        percent_tested,
        test_result_lag,
        sensitivity,
        specificity,
        detection_threshold,
        moving_average_lag,
    )

    @info "Calculating OT characteristics"
    OT_chars = calculate_OutbreakThresholdChars(testarr, incarr)

    return @strdict OT_chars
end

function get_scenario_file(spec)
    dirpath = spec.dirpath
    filecontainer = []

    for file in readdir(dirpath)
        push!(filecontainer, joinpath(dirpath, file))
    end

    if length(filecontainer) != 1
        println("Matched $(length(filecontainer)) files, when should be 1")
    end
    return load(filecontainer...)
end

# end