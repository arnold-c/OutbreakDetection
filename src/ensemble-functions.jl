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
    vaccination_coverage_vec,
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
        vaccination_coverage_vec,
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
            vaccination_coverage,
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
                vaccination_coverage
            ),
            time_p,
            nsims,
        )
    end

    return ensemble_spec_vec
end

function run_ensemble_jump_prob(dict_of_ensemble_params; force = false)
    prog = Progress(length(dict_of_ensemble_params))
    @floop for ensemble_params in dict_of_ensemble_params
        @produce_or_load(
            run_jump_prob,
            ensemble_params,
            "$(ensemble_params[:ensemble_spec].dirpath)";
            filename = "ensemble-solution",
            loadfile = false,
            force = force
        )
        next!(prog)
    end
end

"""
    run_jump_prob(ensemble_param_dict)
"""
function run_jump_prob(ensemble_param_dict)
    @unpack ensemble_spec,
    seed,
    quantile_vec,
    outbreak_spec_dict,
    noise_spec_vec,
    outbreak_detection_spec_vec,
    test_spec_vec = ensemble_param_dict

    @unpack state_parameters, dynamics_parameters, time_parameters, nsims =
        ensemble_spec

    @unpack tstep, tlength, trange = time_parameters

    ensemble_seir_arr = zeros(
        Int64, tlength, size(state_parameters.init_states, 1), nsims
    )
    ensemble_change_arr = zeros(
        Int64, tlength, size(state_parameters.init_states, 1), nsims
    )

    n_transitions = (size(state_parameters.init_states, 1) - 1) * 2 + 1

    ensemble_jump_arr = zeros(Int64, tlength, n_transitions, nsims)
    ensemble_beta_arr = zeros(Float64, tlength)

    for k in 1:nsims
        @views seir = ensemble_seir_arr[:, :, k]
        @views change = ensemble_change_arr[:, :, k]
        @views jump = ensemble_jump_arr[:, :, k]

        run_seed = seed + (k - 1)

        seir_mod!(
            @view(ensemble_seir_arr[:, :, k]),
            @view(ensemble_change_arr[:, :, k]),
            @view(ensemble_jump_arr[:, :, k]),
            ensemble_beta_arr,
            state_parameters.init_states,
            Vector{Float64}(undef, n_transitions),
            dynamics_parameters,
            time_parameters;
            type = "stoch",
            seed = run_seed,
        )
    end

    quantile_param_dict = dict_list(
        @dict(ensemble_spec, ensemble_seir_arr, quantiles = quantile_vec)
    )

    summarize_ensemble_jump_prob(quantile_param_dict)

    for dict in outbreak_spec_dict
        dict[:dirpath] = joinpath(
            ensemble_spec.dirpath, dict[:outbreak_spec].dirpath
        )
        dict[:ensemble_spec] = ensemble_spec
        dict[:ensemble_jump_arr] = ensemble_jump_arr
        dict[:noise_spec_vec] = noise_spec_vec
        dict[:outbreak_detection_spec_vec] = outbreak_detection_spec_vec
        dict[:test_spec_vec] = test_spec_vec
    end

    run_define_outbreaks(outbreak_spec_dict)

    return @strdict ensemble_seir_arr ensemble_spec
end

function summarize_ensemble_jump_prob(dict_of_ensemble_params)
    for ensemble_params in dict_of_ensemble_params
        @produce_or_load(
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
"""
function jump_prob_summary(ensemble_param_dict)
    @unpack ensemble_spec, ensemble_seir_arr, quantiles = ensemble_param_dict
    @unpack state_parameters, dynamics_parameters, time_parameters, nsims =
        ensemble_spec

    @unpack beta_force, annual_births_per_k = dynamics_parameters
    @unpack tstep, tlength, trange = time_parameters
    @unpack init_states, init_state_props = state_parameters

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

    return @strdict ensemble_seir_summary caption quantiles
end

function run_define_outbreaks(dict_of_outbreak_spec_params)
    for outbreak_spec_params in dict_of_outbreak_spec_params
        @produce_or_load(
            define_outbreaks,
            outbreak_spec_params,
            "$(outbreak_spec_params[:dirpath])";
            filename = "ensemble-incidence-array",
            loadfile = false
        )
    end
end

function define_outbreaks(incidence_param_dict)
    @unpack ensemble_spec,
    ensemble_jump_arr,
    outbreak_spec,
    noise_spec_vec,
    outbreak_detection_spec_vec,
    test_spec_vec =
        incidence_param_dict

    ensemble_inc_arr = create_inc_infec_arr(
        ensemble_jump_arr, outbreak_spec
    )

    ensemble_scenarios = create_combinations_vec(
        ScenarioSpecification,
        (
            [ensemble_spec],
            [outbreak_spec],
            noise_spec_vec,
            outbreak_detection_spec_vec,
            test_spec_vec,
        ),
    )

    scenario_param_dict = dict_list(
        @dict(
            scenario_spec = ensemble_scenarios,
            ensemble_jump_arr,
            ensemble_inc_arr
        )
    )

    run_OutbreakThresholdChars_creation(scenario_param_dict)

    return @strdict ensemble_inc_arr
end

function run_OutbreakThresholdChars_creation(
    dict_of_OTchars_params
)
    for OTChars_params in dict_of_OTchars_params
        @produce_or_load(
            OutbreakThresholdChars_creation,
            OTChars_params,
            "$(OTChars_params[:scenario_spec].dirpath)";
            filename = "ensemble-scenario",
            loadfile = false
        )
    end
end

function OutbreakThresholdChars_creation(OT_chars_param_dict)
    @unpack scenario_spec, ensemble_inc_arr, ensemble_jump_arr =
        OT_chars_param_dict
    @unpack noise_specification,
    outbreak_specification,
    outbreak_detection_specification,
    individual_test_specification = scenario_spec

    noise_array = create_poisson_noise_arr(ensemble_inc_arr, noise_specification)

    testarr, posoddsarr = create_testing_arrs(
        ensemble_inc_arr,
        noise_array,
        outbreak_detection_specification,
        individual_test_specification,
    )

    OT_chars = calculate_OutbreakThresholdChars(testarr, ensemble_inc_arr)

    return @strdict OT_chars testarr posoddsarr
end

function get_ensemble_file() end

function get_ensemble_file(
    ensemble_spec::EnsembleSpecification, outbreak_spec::OutbreakSpecification
)
    dirpath = joinpath(ensemble_spec.dirpath, outbreak_spec.dirpath)

    return load(collect_ensemble_file("incidence-array", dirpath)...)
end

function get_ensemble_file(spec::EnsembleSpecification)
    return load(collect_ensemble_file("solution", spec.dirpath)...)
end

function get_ensemble_file(spec::EnsembleSpecification, quantile::Int64)
    return load(collect_ensemble_file("quantiles_$(quantile)", spec.dirpath)...)
end

function get_ensemble_file(spec::ScenarioSpecification)
    return load(collect_ensemble_file("scenario", spec.dirpath)...)
end

function collect_ensemble_file(type, dirpath)
    filecontainer = []
    for f in readdir(dirpath)
        match_ensemble_file!(type, dirpath, filecontainer, f)
    end
    if length(filecontainer) != 1
        println("Matched $(length(filecontainer)) files, when should be 1")
    end
    return filecontainer
end

function match_ensemble_file!(criteria, dirpath, container, file)
    if occursin(criteria, file)
        push!(container, joinpath(dirpath, file))
    end
end

# end
