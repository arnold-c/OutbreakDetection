export evaluate_missing_optimizations

function evaluate_missing_optimizations(
        missing_scenarios::StructVector{OptimizationScenario};
        scheduler = :dynamic, #:serial, :greedy, :static, :dynamic
        save_results = true,
        save_checkpoints = false,
        save_checkpoint_num = 5,
        checkpoint_dir::String = "checkpoint",
        checkpoint_output_filename_base = joinpath(
            checkpoint_dir,
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
        prog = ProgressMeter.Progress(n_missing; desc = "Optimizing scenarios: ", showspeed = true)
    end

    checkpoint_num = 1

    ensemble_groups = group_structvector(missing_scenarios, :ensemble_specification)

    for (ensemble_key, ensemble_scenarios) in ensemble_groups
        verbose && println("Ensemble specification: $(ensemble_key.ensemble_specification.label)")

        ensemble_simulation = generate_single_ensemble(
            ensemble_key.ensemble_specification;
            seed = seed
        )

        outbreak_spec_groups = group_structvector(
            ensemble_scenarios,
            :outbreak_specification
        )

        for (outbreak_spec_key, outbreak_spec_scenarios) in outbreak_spec_groups
            verbose && println("\tOutbreak Specification: $(outbreak_spec_key.outbreak_specification)")

            outbreak_thresholds = calculate_outbreak_thresholds(
                ensemble_simulation,
                outbreak_spec_key.outbreak_specification
            )

            # Validate that all simulations have at least one outbreak
            validation_result = validate_all_simulations_have_outbreaks(
                outbreak_thresholds,
                ensemble_key.ensemble_specification,
                outbreak_spec_key.outbreak_specification
            )

            if Try.iserr(validation_result)
                error(Try.unwrap_err(validation_result))
            end

            noise_trim_groups = group_structvector(
                outbreak_spec_scenarios,
                :noise_level, :noise_type_description
            )

            for (noise_trim_key, noise_trim_scenarios) in noise_trim_groups
                verbose && println("\tNoise level: $(noise_trim_key.noise_level)\n\tNoise type: $(noise_trim_key.noise_type_description)")

                # Track vaccination coverage for this noise scenario
                vaccination_coverage = if noise_trim_key.noise_type_description == :static
                    # For static noise, vaccination coverage is not applicable (set to NaN)
                    NaN
                else
                    # For dynamical noise, optimize to find vaccination coverage
                    optim_res = optimize_dynamic_noise_params_wrapper(
                        ensemble_key.ensemble_specification,
                        ensemble_simulation,
                        noise_trim_key.noise_level,
                        dynamic_noise_optimization_parameters;
                        verbose = verbose_noise_optimization,
                        seed = seed
                    )
                    optim_res.location[1]
                end

                noise_vecs = if noise_trim_key.noise_type_description == :static
                    create_noise_vecs(
                        PoissonNoiseSpecification(noise_trim_key.noise_level),
                        ensemble_key.ensemble_specification,
                        ensemble_simulation,
                        seed = seed
                    )
                else
                    recreate_noise_vecs(
                        ensemble_key.ensemble_specification,
                        vaccination_coverage;
                        seed = seed
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
                            opt_scenarios_vec;
                            scheduler = scheduler
                        ) do opt_scenario_sv
                            @assert length(opt_scenario_sv) == 1
                            optimization_scenario = opt_scenario_sv[1]

                            # Seed the RNG to ensure consistent optimization results across different
                            # noise scenarios, as create_test_positive_vecs modifies the global RNG state
                            Random.seed!(seed)


                            optimization_result = threshold_optimization(
                                optimization_scenario,
                                test_positives_container,
                                ensemble_simulation,
                                outbreak_thresholds,
                                threshold_optimization_parameters,
                                vaccination_coverage,
                            )

                            return optimization_result
                        end

                        BangBang.append!!(all_results, results_batch)

                        verbose && println("checkpoint number $checkpoint_num, result batch length: $(length(results_batch))")
                        if save_results && save_checkpoints && isdir(checkpoint_dir)
                            if checkpoint_num % save_checkpoint_num == 0
                                save_optimization_checkpoint(
                                    StructVector(all_results),
                                    checkpoint_dir,
                                    checkpoint_output_filename_base,
                                    checkpoint_num
                                )
                                verbose && @info "Saved checkpoint $checkpoint_num"
                            end
                            checkpoint_num += 1
                        end

                        verbose && ProgressMeter.update!(prog, length(all_results))
                    end
                end
            end
        end
    end

    return StructVector(all_results)
end
