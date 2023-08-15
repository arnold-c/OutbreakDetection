using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using FLoops
using ProgressMeter

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(funsdir("transmission-functions.jl"))
includet(funsdir("cleaning-functions.jl"))
includet(funsdir("SEIR-model.jl"))

function run_ensemble_jump_prob(params_dict; prog = prog)
    for p in params_dict
        @produce_or_load(
            run_jump_prob,
            p,
            datadir(
                "seasonal-infectivity-import",
                "tau-leaping",
                "N_$(p[:N])",
                "r_$(p[:u₀_prop][:r])",
                "nsims_$(p[:nsims])",
                "births_per_k_$(p[:births_per_k])",
                "beta_force_$(p[:β_force])",
                "tmax_$(p[:tmax])",
                "deltat_$(p[:dt])",
            );
            prefix = "SEIR_tau_sol",
            filename = savename(
                p;
                allowedtypes = (Symbol, Dict, String, Real),
                accesses = [
                    :N, :u₀_prop, :nsims, :tmax, :dt, :births_per_k, :β_force
                ],
                expand = ["u₀_prop"],
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
    @unpack N, u₀_prop, transmission_p, time_p, nsims, dt, β_force,
        births_per_k, seed = param_dict
    @unpack s, e, i, r = u₀_prop
    @unpack R₀, σ, γ = transmission_p
    @unpack tstep, tlength, trange = time_p

    u₀ = convert.(Int64, [s * N, e * N, i * N, r * N, N])
    u0_dict = Dict(zip([:S, :E, :I, :R, :N], u₀))

    μ = births_per_k / (1000 * 365)
    β₀ = calculate_beta(R₀, γ, μ, 1, N)
    ε = (1.06 * μ * (R₀ - 1)) / sqrt(N) # Commuter imports - see p210 Keeling & Rohani

    p = (β₀, β_force, σ, γ, μ, ε, R₀)

    ensemble_seir_arr = zeros(Int64, size(u₀, 1), tlength, nsims)
    ensemble_change_arr = zeros(Int64, size(u₀, 1), tlength, nsims)
    ensemble_jump_arr = zeros(Int64, 9, tlength, nsims)

    @floop for k in 1:nsims
        @views seir = ensemble_seir_arr[:, :, k]
        @views change = ensemble_change_arr[:, :, k]
        @views jump = ensemble_jump_arr[:, :, k]

        seir_mod!(seir, change, jump, u₀, p, trange; dt = tstep, seed = seed)
    end

    return @strdict ensemble_seir_arr ensemble_change_arr ensemble_jump_arr u0_dict param_dict
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
                "r_$(p[:u₀_prop][:r])",
                "nsims_$(p[:nsims])",
                "births_per_k_$(p[:births_per_k])",
                "beta_force_$(p[:β_force])",
                "tmax_$(p[:tmax])",
                "deltat_$(p[:dt])",
            ),
            prefix = "SEIR_tau_quants";
            filename = savename(
                p;
                allowedtypes = (Symbol, Dict, String, Real),
                accesses = [
                    :N, :u₀_prop, :nsims, :tmax, :dt, :births_per_k, :β_force,
                    :quantiles,
                ],
                expand = ["u₀_prop"],
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
    @unpack N, u₀_prop, nsims, dt, tmax, β_force, births_per_k, quantiles =
        param_dict
    @unpack s, e, i, r = u₀_prop

    sim_name = savename(
        "SEIR_tau_sol",
        param_dict,
        "jld2";
        allowedtypes = (Symbol, Dict, String, Real),
        accesses = [:N, :u₀_prop, :nsims, :tmax, :dt, :births_per_k, :β_force],
        expand = ["u₀_prop"],
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
            "beta_force_$β_force",
            "tmax_$tmax",
            "deltat_$dt",
        ),
        sim_name,
    )

    sol_data = load(sim_path)
    @unpack ensemble_seir_arr, u0_dict = sol_data
    S = u0_dict[:S]
    I = u0_dict[:I]
    R = u0_dict[:R]

    qlow = round(0.5 - quantiles / 200; digits = 3)
    qhigh = round(0.5 + quantiles / 200; digits = 3)

    qs = [qlow, 0.5, qhigh]

    ensemble_seir_summary = create_sir_all_sim_quantiles(
        ensemble_seir_arr; quantiles = qs
    )

    caption = "nsims = $nsims, N = $N, S = $S, I = $I, R = $R, β_force = $β_force,\nbirths per k/annum = $births_per_k dt = $dt, quantile int = $quantiles"

    return @strdict ensemble_seir_summary caption u0_dict param_dict
end
