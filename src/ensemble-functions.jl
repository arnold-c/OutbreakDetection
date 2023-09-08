# module EnsembleFunctions
#
# export run_ensemble_jump_prob, run_jump_prob, summarize_ensemble_jump_prob,
#     jump_prob_summary, get_ensemble_file

using DrWatson
using UnPack
using FLoops
using ProgressMeter

# include("transmission-functions.jl")
# # using .TransmissionFunctions
#
# include("cleaning-functions.jl")
# # using .CleaningFunctions
#
# include("SEIR-model.jl")
# # using .SEIRModel
#
# include("structs.jl")
# using .ODStructs

function create_combinations_vec(custom_function, combinations)
    combs = Iterators.product(combinations...)

    return vec(map(combination -> custom_function(combination...), combs))
end

function create_ensemble_spec_combinations(
    beta_force_vec,
    sigma_vec,
    gamma_vec,
    annual_births_per_k_vec,
    R_0_vec,
    N_vec,
    init_states_prop_dict,
    model_types_vec,
    time_p_vec,
    nsims_vec,
)
    ensemble_spec_combinations = Iterators.product(
        beta_force_vec,
        sigma_vec,
        gamma_vec,
        annual_births_per_k_vec,
        R_0_vec,
        N_vec,
        init_states_prop_dict,
        model_types_vec,
        time_p_vec,
        nsims_vec,
    )

    ensemble_spec_vec = Vector(undef, length(ensemble_spec_combinations))

    for (
        i,
        (
            beta_force,
            sigma,
            gamma,
            annual_births_per_k,
            R_0,
            N,
            init_states_prop,
            model_type,
            time_p,
            nsims,
        ),
    ) in enumerate(ensemble_spec_combinations)
        mu = calculate_mu(annual_births_per_k)
        beta_mean = calculate_beta(R_0, gamma, mu, 1, N)
        epsilon = calculate_import_rate(mu, R_0, N)

        ensemble_spec_vec[i] = EnsembleSpecification(
            model_type,
            StateParameters(
                N, init_states_prop
            ),
            DynamicsParameters(
                beta_mean,
                beta_force,
                sigma,
                gamma,
                mu,
                annual_births_per_k,
                epsilon,
                R_0,
            ),
            time_p,
            nsims,
        )
    end

    return ensemble_spec_vec
end

function run_ensemble_jump_prob(dict_of_ensemble_params; prog = prog)
    for ensemble_params in dict_of_ensemble_params
        @produce_or_load(
            run_jump_prob,
            ensemble_params,
            "$(ensemble_params[:ensemble_spec].dirpath)";
            filename = "ensemble-solution",
            loadfile = false
        )
        next!(prog)
    end
end

"""
    run_jump_prob(ensemble_param_dict)
"""
function run_jump_prob(ensemble_param_dict)
    @unpack ensemble_spec, seed = ensemble_param_dict
    @unpack state_parameters, dynamics_parameters, time_parameters, nsims =
        ensemble_spec

    @unpack tstep, tlength, trange = time_parameters

    ensemble_seir_arr = zeros(
        Int64, size(state_parameters.init_states, 1), tlength, nsims
    )
    ensemble_change_arr = zeros(
        Int64, size(state_parameters.init_states, 1), tlength, nsims
    )
    ensemble_jump_arr = zeros(Int64, 9, tlength, nsims)

    @floop for k in 1:nsims
        @views seir = ensemble_seir_arr[:, :, k]
        @views change = ensemble_change_arr[:, :, k]
        @views jump = ensemble_jump_arr[:, :, k]

        run_seed = seed + (k - 1)

        seir_mod!(
            seir,
            change,
            jump,
            state_parameters.init_states,
            dynamics_parameters,
            time_parameters;
            seed = run_seed,
        )
    end

    return @strdict ensemble_seir_arr ensemble_change_arr ensemble_jump_arr dynamics_parameters state_parameters time_parameters ensemble_param_dict
end

function summarize_ensemble_jump_prob(dict_of_ensemble_params; prog = prog)
    for ensemble_params in dict_of_ensemble_params
        @produce_or_load(
            jump_prob_summary,
            ensemble_params,
            "$(ensemble_params[:ensemble_spec].dirpath)";
            filename = "ensemble-quantiles_$(ensemble_params[:quantiles])",
            loadfile = false
        )
        next!(prog)
    end
end

"""
    jump_prob_summary(param_dict)
"""
function jump_prob_summary(ensemble_param_dict)
    @unpack ensemble_spec, seed, quantiles = ensemble_param_dict
    @unpack state_parameters, dynamics_parameters, time_parameters, nsims =
        ensemble_spec

    @unpack beta_force, annual_births_per_k = dynamics_parameters
    @unpack tstep, tlength, trange = time_parameters
    @unpack init_states, init_state_props = state_parameters

    ensemble_sol_file = get_ensemble_file("solution", ensemble_spec)

    @unpack ensemble_seir_arr, state_parameters = ensemble_sol_file
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

    caption = "nsims = $nsims, N = $N, S = $S_init, I = $I_init, R = $R_init, beta_force = $beta_force,\nbirths per k/annum = $annual_births_per_k tstep = $(time_parameters.tstep), quantile int = $quantiles"

    return @strdict ensemble_seir_summary caption state_parameters ensemble_param_dict
end

function get_ensemble_file(type, spec)
    dirpath = spec.dirpath
    filecontainer = []
    for f in readdir(dirpath)
        match_ensemble_file!(type, dirpath, filecontainer, f)
    end
    if length(filecontainer) != 1
        println("Matched $(length(filecontainer)) files, when should be 1")
    end
    if type != "solution"
        try
            parse(Int, type)
        catch
            println("The quantile could not be parsed correctly")
        end
    end
    if type != "solution" && length(filecontainer) == 0
        println(
            "It looks like are trying to return a quantile file. Check that the quantile simulation has been run",
        )
    end
    return load(filecontainer...)
end

function match_ensemble_file!(criteria, dirpath, container, file)
    if occursin(criteria, file)
        push!(container, joinpath(dirpath, file))
    end
end

# end
