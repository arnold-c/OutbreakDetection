export evaluate_missing_optimizations

function evaluate_missing_optimizations(
        missing_scenarios::StructVector{OptimizationScenario};
        scheduler = :dynamic, #:serial, :greedy, :static, :dynamic
        save_results = true,
        save_checkpoints = false,
        save_checkpoint_num = 5,
        checkpoint_dir::String = "",
        checkpoint_output_filename_base = joinpath(
            filedir,
            string(Dates.now()) * "_" * "checkpoint_batch_",
        ),
        verbose::Bool = true,
        verbose_noise_optimization = false,
        seed = 1234,
        dynamic_noise_optimization_parameters::NoiseVaccinationOptimizationParameters = NoiseVaccinationOptimizationParameters(),
        threshold_optimization_parameters::ThresholdOptimizationParameters = ThresholdOptimizationParameters(),
    )
    @assert scheduler in [:dynamic, :static, :greedy, :serial]

    n_missing = length(missing_scenarios)
    n_missing == 0 && return StructVector(OptimizationResult[])

    # Intialize with empty 0-length vector as could be unspecified number of quantiles
    # for a given scenario
    all_results = OptimizationResult[]

    if verbose
        prog = ProgressMeter.Progress(n_missing; desc = "Evaluating grid points: ", showspeed = true)
    end

    checkpoint_num = 1

    ensemble_groups = group_structvector(missing_scenarios, :ensemble_specification)

    for (ensemble_key, ensemble_scenarios) in ensemble_groups
        verbose && println("Ensemble specification: $(ensemble_key.ensemble_specification.label)")

        ensemble_simulation = generate_single_ensemble(
            ensemble_key.ensemble_specification;
            seed = seed
        )
        # TODO: Add grouping based on outbreakdetection specification and create vecs of outbreak statuses
        enddates_vec = fill(
            ensemble_key.ensemble_specification.time_parameters.tlength,
            ensemble_key.ensemble_specification.nsims
        )

        noise_trim_groups = group_structvector(
            ensemble_scenarios,
            :noise_level, :noise_type_description
        )

        for (noise_trim_key, noise_trim_scenarios) in noise_trim_groups
            verbose && println("\tNoise level: $(noise_trim_key.noise_level)\n\tNoise type: $(noise_trim_key.noise_type_description)")

            noise_vecs = if noise_trim_key.noise_type_description == :static
                create_noise_vecs(
                    PoissonNoise(noise_trim_key.noise_level),
                    ensemble_key.ensemble_specification,
                    enddates_vec,
                    ensemble_simulation,
                    seed = seed
                )
            else
                # Call the original optimization function with the computed mean
                optim_res = optimize_dynamic_noise_params_wrapper(
                    ensemble_key.ensemble_specification,
                    ensemble_simulation,
                    enddates_vec,
                    noise_trim_key.noise_level,
                    dynamic_noise_optimization_parameters;
                    verbose = verbose_noise_optimization,
                    seed = seed
                )

                recreate_noise_vecs(
                    ensemble_key.ensemble_specification,
                    enddates_vec,
                    optim_res.location[1],
                )
            end

            test_groups = group_structvector(
                noise_trim_scenarios,
                :percent_tested, :test_specification
            )

            for (test_key, test_scenarios) in test_groups
                verbose && println("\t\tTest Specification: $(test_key.test_specification)\n\t\tTest Percentage: $(test_key.percent_tested)")

                test_positives = create_test_positive_vecs(
                    ensemble_simulation,
                    noise_vecs,
                    test_key.percent_tested,
                    test_key.test_specification,
                )

                alert_method_groups = group_structvector(
                    test_scenarios,
                    :alert_method
                )

                for (alert_method_key, alert_method_scenarios) in alert_method_groups
                    verbose && println("\t\t\tAlert Method: $(alert_method_key.alert_method)")

                    test_positives_container = create_test_positive_container(
                        alert_method_key.alert_method,
                        test_positives
                    )

                    opt_groups = group_structvector(
                        alert_method_scenarios,
                        :accuracy_metric,
                        :threshold_bounds
                    )

                    opt_scenarios_vec = collect(values(opt_groups))

                    results_batch = OhMyThreads.tmap(
                        opt_scenarios_vec; scheduler = scheduler
                    ) do opt_scenario_sv
                        @assert length(opt_scenario_sv) == 1
                        optimization_scenario = opt_scenario_sv[1]

                        optimization_result = threshold_optimization(
                            optimization_scenario,
                            test_positives_container,
                            threshold_optimization_parameters,
                        )

                        return optimization_result
                    end

                    BangBang.append!!(all_results, results_batch)

                    if save_results && save_checkpoints && checkpoint_num % save_checkpoint_num == 0 && !isempty(checkpoint_dir)
                        save_checkpoint_structvector(
                            StructVector(all_results),
                            checkpoint_dir,
                            checkpoint_output_filename_base,
                            checkpoint_num
                        )
                        verbose && @info "Saved checkpoint $checkpoint_num"
                        checkpoint_num += 1
                    end

                    verbose && ProgressMeter.update!(prog, length(all_results))
                end
            end
        end
    end

    return StructVector(all_results)
end
