export run_optimization, optimization_wrapper

"""
    run_optimization(objective_function, OT_chars_param_dict, 
                    optim_method=MSO, accuracy_measure=arithmetic_mean, 
                    kwargs...)

Run threshold optimization for a scenario.

Returns optim_minimizer and optim_minimum.
"""
function run_optimization(
        objective_function,
        OT_chars_param_dict,
        optim_method::OptimizationMethods = MSO,
        accuracy_measure = arithmetic_mean,
        kwargs...,
    )
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

"""
    optimization_wrapper(objective_function_closure, optim_method; kwargs...)

Wrapper for optimization methods.
"""
function optimization_wrapper(
        objective_function_closure,
        optim_method::OptimizationMethods;
        kwargs...,
    )
    return optimization_wrapper(
        objective_function_closure,
        LightSumTypes.variant(optim_method);
        kwargs...
    )
end

"""
    optimization_wrapper(objective_function_closure, optim_method::MSO; kwargs...)

Multistart optimization wrapper using NLopt.
"""
function optimization_wrapper(
        objective_function_closure,
        optim_method::MSO;
        lowers = [0.0],
        uppers = [50.0],
        local_algorithm = NLopt.LN_BOBYQA,
        n_sobol_points = 100,
        use_threads = false,
        xtol_rel = 1.0e-3,
        xtol_abs = 1.0e-3,
        kwargs...,
    )
    kwargs_dict = Dict{Symbol, Any}(kwargs)

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
