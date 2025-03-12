using StatsBase: StatsBase
using NaNMath: NaNMath
using UnPack: @unpack
# using Optim: Optim
using DataFrames: DataFrames
using QuadDIRECT: QuadDIRECT
using MultistartOptimization: MultistartOptimization
using NLopt: NLopt

function setup_optimization(ensemble_param_dict)
    UnPack.@unpack ensemble_spec,
    seed,
    outbreak_spec = ensemble_param_dict

    UnPack.@unpack state_parameters,
    dynamics_parameters, time_parameters,
    nsims =
        ensemble_spec

    UnPack.@unpack tstep, tlength, trange = time_parameters

    ensemble_seir_vecs = Array{typeof(state_parameters.init_states),2}(
        undef,
        tlength,
        nsims,
    )

    ensemble_inc_vecs = Array{typeof(StaticArrays.SVector(0)),2}(
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

    ensemble_inc_arr, ensemble_thresholds_vec = create_inc_infec_arr(
        ensemble_inc_vecs, outbreak_spec
    )

    return ensemble_inc_arr, ensemble_thresholds_vec
end

function run_optimization(
    objective_function,
    OT_chars_param_dict,
    optim_method::TMethod = QD;
    kwargs...,
) where {TMethod<:Type{<:OptimizationMethods}}
    UnPack.@unpack scenario_spec, ensemble_inc_arr, thresholds_vec, seed =
        OT_chars_param_dict
    UnPack.@unpack noise_specification,
    outbreak_detection_specification,
    individual_test_specification = scenario_spec

    noise_array = create_noise_arr(
        noise_specification,
        ensemble_inc_arr;
        ensemble_specification = scenario_spec.ensemble_specification,
        seed = seed,
    )[1]

    obj_inputs = (;
        ensemble_inc_arr,
        noise_array,
        outbreak_detection_specification,
        individual_test_specification,
        thresholds_vec,
    )

    objective_function_closure = x -> objective_function(x, obj_inputs)

    optim_minimizer, optim_minimum = optimization_wrapper(
        objective_function_closure,
        optim_method;
        kwargs...,
    )

    return optim_minimizer, optim_minimum
end

function optimization_wrapper(
    objective_function_closure,
    ::Type{QD};
    splits = ([8.0, 15.0, 35.0],),
    lowers = [0.0],
    uppers = [50.0],
    kwargs...,
)
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    if haskey(kwargs_dict, :splits)
        splits = kwargs_dict[:splits]
    end
    if haskey(kwargs_dict, :lowers)
        lowers = kwargs_dict[:lowers]
    end
    if haskey(kwargs_dict, :uppers)
        uppers = kwargs_dict[:uppers]
    end

    allowed_kwargs = (
        :rtol,
        :atol,
        :fvalue,
        :maxevals,
        :nquasinewton,
        :minwidth,
        :print_interval,
    )

    filtered_kwargs = filter(((k, v),) -> k in allowed_kwargs, kwargs_dict)

    optim_minimizer, optim_minimum = QuadDIRECT.minimize(
        objective_function_closure,
        splits,
        lowers,
        uppers;
        filtered_kwargs...,
    )

    @assert length(optim_minimizer) == 1

    return optim_minimizer[1], optim_minimum
end

function optimization_wrapper(
    objective_function_closure,
    ::Type{MSO};
    lowers = [0.0],
    uppers = [50.0],
    local_algorithm = NLopt.LN_BOBYQA,
    n_sobol_points = 100,
    use_threads = false,
    xtol_rel = 1e-3,
    xtol_abs = 1e-3,
    kwargs...,
)
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    if haskey(kwargs_dict, :lowers)
        lowers = kwargs_dict[:lowers]
    end
    if haskey(kwargs_dict, :uppers)
        uppers = kwargs_dict[:uppers]
    end
    if haskey(kwargs_dict, :n_sobol_points)
        n_sobol_points = kwargs_dict[:n_sobol_points]
    end
    if haskey(kwargs_dict, :use_threads)
        use_threads = kwargs_dict[:use_threads]
    end
    if !haskey(kwargs_dict, :xtol_rel)
        kwargs_dict[:xtol_rel] = xtol_rel
    end
    if !haskey(kwargs_dict, :xtol_abs)
        kwargs_dict[:xtol_abs] = xtol_abs
    end

    allowed_kwargs = (
        :stopval,
        :xtol_rel,
        :xtol_abs,
        :ftol_rel,
        :ftol_abs,
        :maxeval,
        :maxtime,
    )

    filtered_kwargs = filter(((k, v),) -> k in allowed_kwargs, kwargs_dict)

    P = MultistartOptimization.MinimizationProblem(
        objective_function_closure,
        lowers,
        uppers,
    )

    local_method = MultistartOptimization.NLoptLocalMethod(
        local_algorithm;
        filtered_kwargs...,
    )

    multistart_method = MultistartOptimization.TikTak(n_sobol_points)

    p = MultistartOptimization.multistart_minimization(
        multistart_method,
        local_method,
        P;
        use_threads = use_threads,
    )

    @assert length(p.location) == 1

    return p.location[1], p.value
end

function objective_function(
    alert_threshold_vec,
    inputs,
)
    @assert length(alert_threshold_vec) == 1

    @unpack ensemble_inc_arr,
    noise_array,
    outbreak_detection_specification,
    individual_test_specification,
    thresholds_vec = inputs

    outbreak_detection_specification = OutbreakDetectionSpecification(
        alert_threshold_vec[1],
        outbreak_detection_specification.moving_average_lag,
        outbreak_detection_specification.percent_visit_clinic,
        outbreak_detection_specification.percent_clinic_tested,
        outbreak_detection_specification.alert_method.method_name,
    )

    testarr = create_testing_arrs(
        ensemble_inc_arr,
        noise_array,
        outbreak_detection_specification,
        individual_test_specification,
    )[1]

    objective = calculate_ensemble_objective_metric(
        testarr, ensemble_inc_arr, thresholds_vec
    )

    return objective
end

function calculate_ensemble_objective_metric(
    testarr, infecarr, thresholds_vec
)
    mean_accuracy = map(axes(infecarr, 3)) do sim
        dailychars = calculate_daily_detection_characteristics(
            @view(testarr[:, 6, sim]), @view(infecarr[:, 3, sim])
        )
        alertrle = StatsBase.rle(@view(testarr[:, 6, sim]))
        outbreakbounds = thresholds_vec[sim]
        alertbounds = calculate_outbreak_thresholds(alertrle; ncols = 3)

        accuracy = calculate_outbreak_detection_accuracy(
            outbreakbounds, alertbounds
        )
    end

    return 1 - NaNMath.mean(mean_accuracy)
end

function calculate_outbreak_detection_accuracy(
    outbreakbounds, alertbounds
)
    filtered_matched_bounds = match_outbreak_detection_bounds(
        outbreakbounds, alertbounds
    )[1]

    noutbreaks = size(outbreakbounds, 1)
    nalerts = size(alertbounds, 1)

    n_true_outbreaks_detected = length(
        Set(@view(filtered_matched_bounds[:, 1]))
    )
    n_correct_alerts = size(filtered_matched_bounds, 1)

    perc_true_outbreaks_detected = n_true_outbreaks_detected / noutbreaks
    perc_alerts_correct = n_correct_alerts / nalerts # c.f. PPV

    accuracy = NaNMath.mean([perc_true_outbreaks_detected, perc_alerts_correct])

    return accuracy
end
