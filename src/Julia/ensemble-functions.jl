using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using FLoops
using ProgressMeter

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(funsdir("transmission-functions.jl"))
includet(funsdir("cleaning-functions.jl"))
includet(funsdir("SEIR-model.jl"))
includet(funsdir("structs.jl"))

function run_ensemble_jump_prob(params_dict; prog = prog)
    for p in params_dict
        @produce_or_load(
            run_jump_prob,
            p,
            datadir(
                "seasonal-infectivity-import",
                "tau-leaping",
                "N_$(p[:N])",
                "r_$(p[:init_states_prop][:r])",
                "nsims_$(p[:nsims])",
                "births_per_k_$(p[:births_per_k])",
                "beta_force_$(p[:beta_force])",
                "tmax_$(p[:tmax])",
                "tstep_$(p[:tstep])",
            );
            prefix = "SEIR_tau_sol",
            filename = savename(
                p;
                allowedtypes = (Symbol, Dict, String, Real),
                accesses = [
                    :N, :init_states_prop, :nsims, :tmax, :tstep, :births_per_k,
                    :beta_force,
                ],
                expand = ["init_states_prop"],
                sort = false,
            ),
            loadfile = false
        )
        next!(prog)
    end
end

"""
    run_jump_prob(param_dict)
"""
function run_jump_prob(param_dict)
    @unpack N, init_states_prop, transmission_p, time_p, nsims, tstep,
    beta_force,
    births_per_k, seed = param_dict
    @unpack s_init_prop, e_init_prop, i_init_prop, r_init_prop =
        init_states_prop
    @unpack R_0, sigma, gamma = transmission_p
    @unpack tstep, tlength, trange = time_p

    init_states =
        convert.(
            Int64,
            [
                s_init_prop * N,
                e_init_prop * N,
                i_init_prop * N,
                r_init_prop * N,
                N,
            ],
        )
    init_states_dict = Dict(zip([:S, :E, :I, :R, :N], init_states))

    mu = births_per_k / (1_000 * 365)
    beta_mean = calculate_beta(R_0, gamma, mu, 1, N)
    epsilon = (1.06 * mu * (R_0 - 1)) / sqrt(N) # Commuter imports - see p210 Keeling & Rohani

    p = (beta_mean, beta_force, sigma, gamma, mu, epsilon, R_0)

    ensemble_seir_arr = zeros(Int64, size(init_states, 1), tlength, nsims)
    ensemble_change_arr = zeros(Int64, size(init_states, 1), tlength, nsims)
    ensemble_jump_arr = zeros(Int64, 9, tlength, nsims)

    @floop for k in 1:nsims
        @views seir = ensemble_seir_arr[:, :, k]
        @views change = ensemble_change_arr[:, :, k]
        @views jump = ensemble_jump_arr[:, :, k]

        seir_mod!(
            seir,
            change,
            jump,
            init_states,
            p,
            trange;
            tstep = tstep,
            seed = seed,
        )
    end

    return @strdict ensemble_seir_arr ensemble_change_arr ensemble_jump_arr init_states_dict param_dict
end

function summarize_ensemble_jump_prob(params_dict; prog = prog)
    for p in params_dict
        @produce_or_load(
            jump_prob_summary,
            p,
            datadir(
                "seasonal-infectivity-import",
                "tau-leaping",
                "N_$(p[:N])",
                "r_$(p[:init_states_prop][:r])",
                "nsims_$(p[:nsims])",
                "births_per_k_$(p[:births_per_k])",
                "beta_force_$(p[:beta_force])",
                "tmax_$(p[:tmax])",
                "tstep_$(p[:tstep])",
            ),
            prefix = "SEIR_tau_quants";
            filename = savename(
                p;
                allowedtypes = (Symbol, Dict, String, Real),
                accesses = [
                    :N, :init_states_prop, :nsims, :tmax, :tstep, :births_per_k,
                    :beta_force,
                    :quantiles,
                ],
                expand = ["init_states_prop"],
                sort = false,
            ),
            loadfile = false
        )
        next!(prog)
    end
end

"""
    jump_prob_summary(param_dict)
"""
function jump_prob_summary(param_dict)
    @unpack N,
    init_states_prop, nsims, tstep, tmax, beta_force, births_per_k,
    quantiles =
        param_dict
    @unpack s_init_prop, e_init_prop, i_init_prop, r_init_prop =
        init_states_prop

    sim_name = savename(
        "SEIR_tau_sol",
        param_dict,
        "jld2";
        allowedtypes = (Symbol, Dict, String, Real),
        accesses = [
            :N,
            :init_states_prop,
            :nsims,
            :tmax,
            :tstep,
            :births_per_k,
            :beta_force,
        ],
        expand = ["init_states_prop"],
        sort = false,
    )
    sim_path = joinpath(
        datadir(
            "seasonal-infectivity-import",
            "tau-leaping",
            "N_$N",
            "r_$r",
            "nsims_$nsims",
            "births_per_k_$births_per_k",
            "beta_force_$beta_force",
            "tmax_$tmax",
            "tstep_$dt",
        ),
        sim_name,
    )

    sol_data = load(sim_path)
    @unpack ensemble_seir_arr, init_states_dict = sol_data
    S_init = init_states_dict[:S]
    I_init = init_states_dict[:I]
    R_init = init_states_dict[:R]

    qlow = round(0.5 - quantiles / 200; digits = 3)
    qhigh = round(0.5 + quantiles / 200; digits = 3)

    qs = [qlow, 0.5, qhigh]

    ensemble_seir_summary = create_sir_all_sim_quantiles(
        ensemble_seir_arr; quantiles = qs
    )

    caption = "nsims = $nsims, N = $N, S = $S_init, I = $I_init, R = $R_init, beta_force = $beta_force,\nbirths per k/annum = $births_per_k tstep = $tstep, quantile int = $quantiles"

    return @strdict ensemble_seir_summary caption init_states_dict param_dict
end

function get_ensemble_file(type, spec)
    dirpath = get_ensemble_file_dir(spec)
    filecontainer = []
    for f in readdir(dirpath)
        match_ensemble_file!(type, dirpath, filecontainer, f)
    end
    if length(filecontainer) != 1
        println("Matched $(length(filecontainer)) files, when should be 1")
    end
    if type != "sol"
        try
            parse(Int, type)
        catch
            println("The quantile could not be parsed correctly")
        end
    end
    if type != "sol" && length(filecontainer) == 0
        println(
            "It looks like are trying to return a quantile file. Check that the quantile simulation has been run",
        )
    end
    return load(filecontainer...)
end

function get_ensemble_file_dir(spec)
    dirpath = joinpath(
        spec.modeltypes...,
        "N_$(spec.N)",
        "r_$(spec.Rinit_prop)",
        "nsims_$(spec.nsims)",
        "births_per_k_$(spec.births_per_k)",
        "beta_force_$(spec.beta_force)",
        "tmax_$(spec.time_parameters.tmax)",
        "tstep_$(spec.time_parameters.tstep)",
    )
    return datadir(dirpath)
end

function match_ensemble_file!(criteria, dirpath, container, file)
    if occursin(criteria, file)
        push!(container, joinpath(dirpath, file))
    end
end