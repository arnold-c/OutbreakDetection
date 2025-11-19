using DrWatson: DrWatson
using UnPack: UnPack
using FLoops: FLoops
using ProgressMeter: ProgressMeter
using StaticArrays: StaticArrays

export run_ensemble_jump_prob, run_jump_prob,
    summarize_ensemble_jump_prob, jump_prob_summary

"""
    run_ensemble_jump_prob(dict_of_ensemble_params; force=false)

Run ensemble simulations for multiple parameter sets.
"""
function run_ensemble_jump_prob(dict_of_ensemble_params; force = false)
    prog = ProgressMeter.Progress(length(dict_of_ensemble_params))
    for ensemble_params in dict_of_ensemble_params
        DrWatson.@produce_or_load(
            run_jump_prob,
            ensemble_params,
            "$(ensemble_params[:ensemble_spec].dirpath)";
            filename = "ensemble-solution",
            loadfile = false,
            force = force
        )
        ProgressMeter.next!(prog)
    end
    return
end

"""
    run_jump_prob(ensemble_param_dict)

Run a single ensemble simulation.
"""
function run_jump_prob(ensemble_param_dict)
    UnPack.@unpack ensemble_spec,
        seed,
        quantile_vec,
        outbreak_spec_dict,
        noise_spec_vec,
        outbreak_detection_spec_vec,
        test_spec_vec = ensemble_param_dict

    UnPack.@unpack state_parameters,
        dynamics_parameters, time_parameters,
        nsims =
        ensemble_spec

    UnPack.@unpack tstep, tlength, trange = time_parameters

    ensemble_seir_vecs = Array{typeof(state_parameters.init_states), 2}(
        undef,
        tlength,
        nsims,
    )

    ensemble_inc_vecs = Array{typeof(StaticArrays.SVector(0)), 2}(
        undef,
        tlength,
        nsims,
    )

    ensemble_beta_arr = zeros(Float64, tlength)

    for sim in axes(ensemble_inc_vecs, 2)
        run_seed = seed + (sim - 1)

        seir_mod!(
            @view(ensemble_seir_vecs[:, sim]),
            @view(ensemble_inc_vecs[:, sim]),
            ensemble_beta_arr,
            state_parameters.init_states,
            dynamics_parameters,
            time_parameters;
            seed = run_seed,
        )
    end

    ensemble_seir_arr = convert_svec_to_array(ensemble_seir_vecs)

    quantile_param_dict = DrWatson.dict_list(
        DrWatson.@dict(
            ensemble_spec, ensemble_seir_arr, quantiles = quantile_vec
        )
    )

    summarize_ensemble_jump_prob(quantile_param_dict)

    for dict in outbreak_spec_dict
        dict[:dirpath] = joinpath(
            ensemble_spec.dirpath, dict[:outbreak_spec].dirpath
        )
        dict[:ensemble_spec] = ensemble_spec
        dict[:ensemble_inc_vecs] = ensemble_inc_vecs
        dict[:noise_spec_vec] = noise_spec_vec
        dict[:outbreak_detection_spec_vec] = outbreak_detection_spec_vec
        dict[:test_spec_vec] = test_spec_vec
        dict[:seed] = seed
    end

    run_define_outbreaks(outbreak_spec_dict)

    return DrWatson.@strdict ensemble_seir_arr ensemble_spec
end

"""
    summarize_ensemble_jump_prob(dict_of_ensemble_params)

Summarize ensemble simulation results with quantiles.
"""
function summarize_ensemble_jump_prob(dict_of_ensemble_params)
    return FLoops.@floop for ensemble_params in dict_of_ensemble_params
        DrWatson.@produce_or_load(
            jump_prob_summary,
            ensemble_params,
            "$(ensemble_params[:ensemble_spec].dirpath)";
            filename = "ensemble-quantiles_$(ensemble_params[:quantiles])",
            loadfile = false
        )
    end
end

"""
    jump_prob_summary(param_dict)

Create summary statistics for ensemble simulation.
"""
function jump_prob_summary(ensemble_param_dict)
    UnPack.@unpack ensemble_spec, ensemble_seir_arr, quantiles =
        ensemble_param_dict
    UnPack.@unpack state_parameters,
        dynamics_parameters, time_parameters,
        nsims =
        ensemble_spec

    UnPack.@unpack beta_force, annual_births_per_k = dynamics_parameters
    UnPack.@unpack tstep, tlength, trange = time_parameters
    UnPack.@unpack init_states, init_state_props = state_parameters

    N = init_states[:N]
    S_init = init_states[:S]
    I_init = init_states[:I]
    R_init = init_states[:R]

    qlow = round(0.5 - quantiles / 200; digits = 3)
    qhigh = round(0.5 + quantiles / 200; digits = 3)

    qs = [qlow, 0.5, qhigh]

    ensemble_seir_summary = create_sir_all_sim_quantiles(
        ensemble_seir_arr; quantiles = qs
    )

    caption = "nsims = $nsims, N = $N, S = $S_init, I = $I_init, R = $R_init, beta_force = $beta_force,\nbirths per k/annum = $annual_births_per_k, tstep = $(time_parameters.tstep), quantile int = $quantiles"

    return DrWatson.@strdict ensemble_seir_summary caption quantiles
end
