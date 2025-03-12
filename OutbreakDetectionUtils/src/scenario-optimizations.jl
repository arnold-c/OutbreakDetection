using DataFrames: DataFrames
using DrWatson: @dict
using StructArrays: StructVector
using ProgressMeter: Progress, next!
using FLoops: FLoops

function run_scenario_optimizations(
    ensemble_specifications,
    outbreak_specifications,
    noise_specifications,
    outbreak_detection_specifications,
    individual_test_specifications,
    optim_method::TMethod = QD;
    seed = 1234,
    executor = FLoops.SequentialEx(),
    kwargs...,
) where {TMethod<:Type{<:OptimizationMethods}}
    optim_df = DataFrames.DataFrame((
        ensemble_spec = EnsembleSpecification[],
        outbreak_spec = OutbreakSpecification[],
        noise_spec = NoiseSpecification[],
        outbreak_detection_spec = OutbreakDetectionSpecification[],
        test_sensitivity = Float64[],
        test_specificity = Float64[],
        test_result_lag = Int64[],
        optimal_threshold = Float64[],
        optimal_accuracy = Float64[],
        optimization_method = Union{Type{QD},Type{MSO}}[],
        OT_chars = StructVector{<:OutbreakThresholdChars}[],
    ))

    run_scenario_optimizations!(
        optim_df,
        ensemble_specifications,
        outbreak_specifications,
        noise_specifications,
        outbreak_detection_specifications,
        individual_test_specifications,
        optim_method;
        seed = seed,
        kwargs...,
    )
    return optim_df
end

function run_scenario_optimizations!(
    optim_df,
    ensemble_specifications,
    outbreak_specifications,
    noise_specifications,
    outbreak_detection_specifications,
    individual_test_specifications,
    optim_method::TMethod = QD;
    seed = 1234,
    executor = FLoops.SequentialEx(),
    kwargs...,
) where {TMethod<:Type{<:OptimizationMethods}}
	ncombinations = length(
            Iterators.product(
                ensemble_specifications,
                outbreak_specifications,
                noise_specifications,
                outbreak_detection_specifications,
                individual_test_specifications,
            ),
        )

    prog = Progress(ncombinations)

	ncompleted = 0
    for ensemble_spec in ensemble_specifications
        for outbreak_spec in outbreak_specifications
            base_param_dict = @dict(
                ensemble_spec = ensemble_spec,
                outbreak_spec = outbreak_spec,
                seed = seed,
            )

            ensemble_inc_arr, thresholds_vec = setup_optimization(
                base_param_dict
            )

            FLoops.@floop executor for noise_spec in noise_specifications
                noise_array, noise_means = create_noise_arr(
                    noise_spec,
                    ensemble_inc_arr;
                    ensemble_specification = ensemble_spec,
                    seed = seed,
                )

                for outbreak_detection_spec in
                    outbreak_detection_specifications,
                    individual_test_spec in individual_test_specifications

					println("Individual Test: $(individual_test_spec.sensitivity), $(individual_test_spec.specificity), $(individual_test_spec.test_result_lag)")

                    obj_inputs = (;
                        ensemble_inc_arr,
                        noise_array,
                        outbreak_detection_specification = outbreak_detection_spec,
                        individual_test_specification = individual_test_spec,
                        thresholds_vec,
                    )

                    objective_function_closure =
                        x -> objective_function(x, obj_inputs)

                    optim_minimizer, optim_minimum = optimization_wrapper(
                        objective_function_closure,
                        optim_method;
                        kwargs...,
                    )

                    optimal_outbreak_detection_spec = OutbreakDetectionSpecification(
                        optim_minimizer,
                        outbreak_detection_spec.moving_average_lag,
                        outbreak_detection_spec.percent_visit_clinic,
                        outbreak_detection_spec.percent_clinic_tested,
                        outbreak_detection_spec.alert_method.method_name,
                    )

                    testarr, test_movingvg_arr = create_testing_arrs(
                        ensemble_inc_arr,
                        noise_array,
                        optimal_outbreak_detection_spec,
                        individual_test_spec,
                    )

                    OT_chars = calculate_OutbreakThresholdChars(
                        testarr, ensemble_inc_arr, thresholds_vec, noise_means
                    )

                    push!(
                        optim_df,
                        (
                            ensemble_spec,
                            outbreak_spec,
                            noise_spec,
                            outbreak_detection_spec,
                            individual_test_spec.sensitivity,
                            individual_test_spec.specificity,
                            individual_test_spec.test_result_lag,
                            optim_minimizer,
                            1 - optim_minimum,
                            optim_method,
                            OT_chars,
                        ),
                    )

					ncompleted += 1
					println("Completed $ncompleted scenarios. $(ncombinations - ncompleted) left of $ncombinations")
                    next!(prog)
					println("")
                end
            end
        end
    end
end
