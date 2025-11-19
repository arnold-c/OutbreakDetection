export calculate_combination_to_run

"""
    calculate_combination_to_run(ensemble_specifications, outbreak_specifications, 
                                 noise_specifications, outbreak_detection_specifications, 
                                 individual_test_specifications, optim_method, 
                                 accuracy_functions; scenario_parameter_symbols=...)

Calculate all combinations of scenarios to run for optimization.
"""
function calculate_combination_to_run(
        ensemble_specifications,
        outbreak_specifications,
        noise_specifications,
        outbreak_detection_specifications,
        individual_test_specifications,
        optim_method::TMethod = MSO,
        accuracy_functions = [arithmetic_mean, calculate_f_beta_score];
        scenario_parameter_symbols = [
            :ensemble_spec,
            :outbreak_spec,
            :noise_spec,
            :outbreak_detection_spec,
            :test_spec,
            :optimization_method,
            :accuracy_function,
        ],
    ) where {TMethod <: Type{<:OptimizationMethods}}
    @assert mapreduce(
        f -> in(f, [arithmetic_mean, calculate_f_beta_score]),
        +,
        unique(accuracy_functions),
    ) == length(unique(accuracy_functions))

    return DataFrames.DataFrame(
        Iterators.product(
            ensemble_specifications,
            outbreak_specifications,
            noise_specifications,
            outbreak_detection_specifications,
            individual_test_specifications,
            [optim_method],
            accuracy_functions,
        ),
        scenario_parameter_symbols,
    )
end
