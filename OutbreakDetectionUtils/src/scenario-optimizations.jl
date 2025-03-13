using DataFrames: DataFrames
using DrWatson: @dict
using StructArrays: StructVector
using ProgressMeter: Progress, next!
using FLoops: FLoops
using Try: Try
using TryExperimental: TryExperimental
using REPL: REPL
using REPL.TerminalMenus: RadioMenu, request
using StyledStrings

function run_scenario_optimizations(
    ensemble_specifications,
    outbreak_specifications,
    noise_specifications,
    outbreak_detection_specifications,
    individual_test_specifications,
    optim_method::TMethod = MSO;
    seed = 1234,
    executor = FLoops.SequentialEx(),
    kwargs...,
) where {TMethod<:Type{<:OptimizationMethods}}
    optim_df = DataFrames.DataFrame((
        ensemble_spec = typeof(ensemble_specifications[1])[],
        outbreak_spec = typeof(outbreak_specifications[1])[],
        noise_spec = Union{
            PoissonNoiseSpecification,DynamicalNoiseSpecification
        }[],
        outbreak_detection_spec = typeof(outbreak_detection_specifications[1])[],
        test_spec = typeof(individual_test_specifications[1])[],
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
    optim_method::TMethod = MSO;
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
                println(
                    styled"{green:\n=================================================================}"
                )
                println(
                    styled"Noise type: {green,inverse: $(getdirpath(noise_spec))}"
                )

                noise_array, noise_means = create_noise_arr(
                    noise_spec,
                    ensemble_inc_arr;
                    ensemble_specification = ensemble_spec,
                    seed = seed,
                )

                for outbreak_detection_spec in
                    outbreak_detection_specifications,
                    individual_test_spec in individual_test_specifications

                    println(
                        styled"\t\t\t\t-> Test specification: {blue: $(get_test_description(individual_test_spec))}, Percent tested: {red,inverse: $(outbreak_detection_spec.percent_tested)}"
                    )

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
                            individual_test_spec,
                            optim_minimizer,
                            1 - optim_minimum,
                            optim_method,
                            OT_chars,
                        ),
                    )

                    next!(prog)
                    println("")
                end
            end
        end
    end
end

function check_missing_scenario_optimizations(
    optim_df,
    ensemble_specifications,
    outbreak_specifications,
    noise_specifications,
    outbreak_detection_specifications,
    individual_test_specifications,
    optim_method::TMethod = MSO;
    disable_time_check = false,
    time_per_run_s = 45,
) where {TMethod<:Type{<:OptimizationMethods}}
    scenario_parameter_symbols = [
        :ensemble_spec,
        :outbreak_spec,
        :noise_spec,
        :outbreak_detection_spec,
        :test_spec,
        :optimization_method,
    ]

    combinations_to_run = DataFrames.DataFrame(
        Iterators.product(
            ensemble_specifications,
            outbreak_specifications,
            noise_specifications,
            outbreak_detection_specifications,
            individual_test_specifications,
            [optim_method],
        ),
        scenario_parameter_symbols,
    )

    missing_combinations = DataFrames.antijoin(
        combinations_to_run,
        optim_df;
        on = scenario_parameter_symbols,
    )

    missing_runs = DataFrames.nrow(missing_combinations)
    if missing_runs == 0
        return Try.Err("No missing simulations")
    end

    if missing_runs > 0 && !disable_time_check
        nrun_time_s = missing_runs * time_per_run_s
        nrun_time_minutes = round(nrun_time_s / 60; digits = 2)
        nrun_time_message = if nrun_time_s < 10
            "less than 10 seconds"
        elseif nrun_time_s < 60
            "approximately $(round(nrun_time_s; digits = 0)) seconds"
        else
            "approximately $(nrun_time_minutes) minutes"
        end
        choice = request(
            "There are $(missing_runs) missing simulations. This is estimated to take $(nrun_time_message). Do you want to continue?",
            RadioMenu(["No", "Yes"]; ctrl_c_interrupt = false),
        )

        if choice != 2
            return Try.Err("User aborted")
        end

        println("Continuing ...")
    end

    return Try.Ok(missing_combinations)
end

function run_missing_scenario_optimizations!(
    optim_df,
    missing_optimizations;
    seed = 1234,
    kwargs...,
)
    missing_optims_df = Try.unwrap(missing_optimizations)

    prog = Progress(DataFrames.nrow(missing_optims_df))

    incidence_sim_grouped_df = DataFrames.groupby(
        missing_optims_df,
        [:ensemble_spec, :outbreak_spec],
    )

    for inc_gp in incidence_sim_grouped_df
        ensemble_spec = inc_gp[1, :ensemble_spec]
        outbreak_spec = inc_gp[1, :outbreak_spec]

        base_param_dict = @dict(
            ensemble_spec = ensemble_spec,
            outbreak_spec = outbreak_spec,
            seed = seed,
        )

        ensemble_inc_arr, thresholds_vec = setup_optimization(
            base_param_dict
        )

        for noise_gp in DataFrames.groupby(inc_gp, [:noise_spec])
            noise_spec = noise_gp[1, :noise_spec]

            noise_array, noise_means = create_noise_arr(
                noise_spec,
                ensemble_inc_arr;
                ensemble_specification = ensemble_spec,
                seed = seed,
            )

            for detect_test_gp in DataFrames.groupby(
                noise_gp, [:outbreak_detection_spec, :test_spec]
            )
                outbreak_detection_spec = detect_test_gp[
                    1, :outbreak_detection_spec
                ]
                individual_test_spec = detect_test_gp[1, :test_spec]

                println(
                    "Individual Test: $(individual_test_spec.sensitivity), $(individual_test_spec.specificity), $(individual_test_spec.test_result_lag)"
                )

                obj_inputs = (;
                    ensemble_inc_arr,
                    noise_array,
                    outbreak_detection_specification = outbreak_detection_spec,
                    individual_test_specification = individual_test_spec,
                    thresholds_vec,
                )

                objective_function_closure =
                    x -> objective_function(x, obj_inputs)

                @assert length(
                    unique(detect_test_gp[:, :optimization_method])
                ) == 1

                optim_method = detect_test_gp[1, :optimization_method]

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

                testarr = create_testing_arrs(
                    ensemble_inc_arr,
                    noise_array,
                    optimal_outbreak_detection_spec,
                    individual_test_spec,
                )[1]

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
                        individual_test_spec,
                        optim_minimizer,
                        1 - optim_minimum,
                        optim_method,
                        OT_chars,
                    ),
                )

                next!(prog)
                println("")
            end
        end
    end

    return nothing
end
