@testset "missing-results.jl" begin
    using OutbreakDetectionCore
    using StructArrays
    using Test
    using Dates

    # Helper function to create minimal test data
    function create_test_ensemble_spec(label::String)

        state_params = StateParameters(;
            N = 100_000,
            s_prop = 0.05,
            e_prop = 0.0,
            i_prop = 0.0
        )

        time_params = SimTimeParameters(;
            tmin = 0.0,
            tmax = 365.0,
            tstep = 1.0
        )

        target_dynamics = TargetDiseaseDynamicsParameters(;
            R_0 = 16.0,
            latent_period = Day(10),
            infectious_duration = Day(8),
            beta_force = 0.2
        )

        common_dynamics = CommonDiseaseDynamicsParameters(;
            births_per_k_pop = 12.0,
            nsims = 10
        )

        dynamics_params = DynamicsParameterSpecification(
            state_params,
            target_dynamics,
            common_dynamics
        )

        noise_params = DynamicalNoiseParameters(;
            R_0 = 0.1,
            latent_period = Day(1),
            infectious_duration = Day(1),
            correlation = "none",
            poisson_component = 0.0
        )

        return EnsembleSpecification(
            label,
            state_params,
            time_params,
            dynamics_params,
            noise_params,
            10  # nsims
        )
    end

    function create_test_scenario(;
            label = "test",
            noise_level = 0.5,
            noise_type = :static,
            percent_tested = 0.5,
            threshold_lower = 1.0,
            threshold_upper = 10.0
        )
        ensemble_spec = create_test_ensemble_spec(label)

        test_spec = IndividualTestSpecification(
            sensitivity = 1.0,
            specificity = 0.8,
            test_result_lag = 0
        )

        alert_method = AlertMethod(DailyThreshold())
        accuracy_metric = AccuracyMetric(F1())

        outbreak_spec = OutbreakSpecification(
            1,  # outbreak_threshold (must be Integer)
            7,  # minimum_outbreak_duration
            10  # minimum_outbreak_size
        )

        return OptimizationScenario(
            ensemble_specification = ensemble_spec,
            noise_level = noise_level,
            noise_type_description = noise_type,
            test_specification = test_spec,
            percent_tested = percent_tested,
            alert_method = alert_method,
            accuracy_metric = accuracy_metric,
            threshold_bounds = (lower = threshold_lower, upper = threshold_upper),
            outbreak_specification = outbreak_spec
        )
    end

    function create_test_result(scenario::OptimizationScenario; optimal_threshold = 5.0)
        return OptimizationResult(
            ensemble_specification = scenario.ensemble_specification,
            noise_level = scenario.noise_level,
            noise_type_description = scenario.noise_type_description,
            test_specification = scenario.test_specification,
            percent_tested = scenario.percent_tested,
            alert_method = scenario.alert_method,
            accuracy_metric = scenario.accuracy_metric,
            threshold_bounds = scenario.threshold_bounds,
            outbreak_specification = scenario.outbreak_specification,
            optimal_threshold = optimal_threshold,
            accuracies = [0.9, 0.85, 0.88],
            proportion_alerts_correct = [0.92, 0.88, 0.9],
            proportion_outbreaks_detected = [0.88, 0.82, 0.86],
            detection_delays = [[5], [6], [4]],
            unavoidable_cases = [[100], [120], [110]],
            alert_durations = [[10], [12], [11]],
            outbreak_durations = [[20], [22], [21]],
            proportion_timeseries_in_alert = [0.1, 0.12, 0.11],
            proportion_timeseries_in_outbreak = [0.2, 0.22, 0.21]
        )
    end

    @testset "scenario_key" begin
        @testset "generates consistent keys for same scenario" begin
            scenario1 = create_test_scenario(label = "test1")
            scenario2 = create_test_scenario(label = "test1")

            key1 = OutbreakDetectionCore.scenario_key(scenario1)
            key2 = OutbreakDetectionCore.scenario_key(scenario2)

            @test key1 == key2
            @test key1 isa Tuple
            @test length(key1) == 9
        end

        @testset "generates different keys for different scenarios" begin
            scenario1 = create_test_scenario(label = "test1", noise_level = 0.5)
            scenario2 = create_test_scenario(label = "test2", noise_level = 0.7)

            key1 = OutbreakDetectionCore.scenario_key(scenario1)
            key2 = OutbreakDetectionCore.scenario_key(scenario2)

            @test key1 != key2
        end

        @testset "key contains all scenario parameters" begin
            scenario = create_test_scenario()
            key = OutbreakDetectionCore.scenario_key(scenario)

            @test key[1] == scenario.ensemble_specification
            @test key[2] == scenario.noise_level
            @test key[3] == scenario.noise_type_description
            @test key[4] == scenario.test_specification
            @test key[5] == scenario.percent_tested
            @test key[6] == scenario.alert_method
            @test key[7] == scenario.accuracy_metric
            @test key[8] == scenario.threshold_bounds
            @test key[9] == scenario.outbreak_specification
        end
    end

    @testset "result_key" begin
        @testset "generates consistent keys for same result" begin
            scenario = create_test_scenario()
            result1 = create_test_result(scenario)
            result2 = create_test_result(scenario)

            key1 = OutbreakDetectionCore.result_key(result1)
            key2 = OutbreakDetectionCore.result_key(result2)

            @test key1 == key2
            @test key1 isa Tuple
            @test length(key1) == 9
        end

        @testset "generates different keys for different results" begin
            scenario1 = create_test_scenario(label = "test1", noise_level = 0.5)
            scenario2 = create_test_scenario(label = "test2", noise_level = 0.7)

            result1 = create_test_result(scenario1)
            result2 = create_test_result(scenario2)

            key1 = OutbreakDetectionCore.result_key(result1)
            key2 = OutbreakDetectionCore.result_key(result2)

            @test key1 != key2
        end

        @testset "result key matches scenario key" begin
            scenario = create_test_scenario()
            result = create_test_result(scenario)

            scenario_k = OutbreakDetectionCore.scenario_key(scenario)
            result_k = OutbreakDetectionCore.result_key(result)

            @test scenario_k == result_k
        end
    end

    @testset "build_scenario_dict" begin
        @testset "creates empty dict for empty results" begin
            empty_results = StructVector{OptimizationResult}(undef, 0)
            dict = OutbreakDetectionCore.build_scenario_dict(empty_results)

            @test dict isa Dict{Tuple, Bool}
            @test isempty(dict)
        end

        @testset "creates dict with correct keys" begin
            scenario1 = create_test_scenario(label = "test1", noise_level = 0.5)
            scenario2 = create_test_scenario(label = "test2", noise_level = 0.7)

            result1 = create_test_result(scenario1)
            result2 = create_test_result(scenario2)

            results = StructVector([result1, result2])
            dict = OutbreakDetectionCore.build_scenario_dict(results)

            @test length(dict) == 2
            @test haskey(dict, OutbreakDetectionCore.result_key(result1))
            @test haskey(dict, OutbreakDetectionCore.result_key(result2))
            @test dict[OutbreakDetectionCore.result_key(result1)] == true
            @test dict[OutbreakDetectionCore.result_key(result2)] == true
        end

        @testset "handles duplicate results" begin
            scenario = create_test_scenario()
            result1 = create_test_result(scenario, optimal_threshold = 5.0)
            result2 = create_test_result(scenario, optimal_threshold = 6.0)

            results = StructVector([result1, result2])
            dict = OutbreakDetectionCore.build_scenario_dict(results)

            # Should only have one entry since keys are the same
            @test length(dict) == 1
            @test haskey(dict, OutbreakDetectionCore.result_key(result1))
        end
    end

    @testset "find_missing_scenarios" begin
        @testset "returns all scenarios when no completed results" begin
            scenario1 = create_test_scenario(label = "test1", noise_level = 0.5)
            scenario2 = create_test_scenario(label = "test2", noise_level = 0.7)
            scenario3 = create_test_scenario(label = "test3", noise_level = 0.9)

            all_scenarios = StructVector([scenario1, scenario2, scenario3])
            empty_results = StructVector{OptimizationResult}(undef, 0)

            missing = find_missing_scenarios(all_scenarios, empty_results)

            @test length(missing) == 3
            @test missing == all_scenarios
        end

        @testset "returns empty when all scenarios completed" begin
            scenario1 = create_test_scenario(label = "test1", noise_level = 0.5)
            scenario2 = create_test_scenario(label = "test2", noise_level = 0.7)

            all_scenarios = StructVector([scenario1, scenario2])

            result1 = create_test_result(scenario1)
            result2 = create_test_result(scenario2)
            completed_results = StructVector([result1, result2])

            missing = find_missing_scenarios(all_scenarios, completed_results)

            @test length(missing) == 0
            @test isempty(missing)
        end

        @testset "identifies partially completed scenarios" begin
            scenario1 = create_test_scenario(label = "test1", noise_level = 0.5)
            scenario2 = create_test_scenario(label = "test2", noise_level = 0.7)
            scenario3 = create_test_scenario(label = "test3", noise_level = 0.9)

            all_scenarios = StructVector([scenario1, scenario2, scenario3])

            # Only complete scenario1 and scenario3
            result1 = create_test_result(scenario1)
            result3 = create_test_result(scenario3)
            completed_results = StructVector([result1, result3])

            missing = find_missing_scenarios(all_scenarios, completed_results)

            @test length(missing) == 1
            @test missing[1].noise_level == 0.7  # scenario2
        end

        @testset "handles single missing scenario" begin
            scenario1 = create_test_scenario(label = "test1", noise_level = 0.5)
            scenario2 = create_test_scenario(label = "test2", noise_level = 0.7)

            all_scenarios = StructVector([scenario1, scenario2])

            result1 = create_test_result(scenario1)
            completed_results = StructVector([result1])

            missing = find_missing_scenarios(all_scenarios, completed_results)

            @test length(missing) == 1
            @test missing[1].noise_level == 0.7
        end

        @testset "handles scenarios with different parameters" begin
            # Create scenarios varying multiple parameters
            scenario1 = create_test_scenario(
                label = "test1",
                noise_level = 0.5,
                percent_tested = 0.5,
                threshold_lower = 1.0
            )
            scenario2 = create_test_scenario(
                label = "test2",
                noise_level = 0.5,
                percent_tested = 0.7,
                threshold_lower = 1.0
            )
            scenario3 = create_test_scenario(
                label = "test3",
                noise_level = 0.7,
                percent_tested = 0.5,
                threshold_lower = 2.0
            )

            all_scenarios = StructVector([scenario1, scenario2, scenario3])

            # Complete only scenario2
            result2 = create_test_result(scenario2)
            completed_results = StructVector([result2])

            missing = find_missing_scenarios(all_scenarios, completed_results)

            @test length(missing) == 2
            @test scenario1 in missing
            @test scenario3 in missing
            @test !(scenario2 in missing)
        end

        @testset "preserves order of missing scenarios" begin
            scenarios = [
                create_test_scenario(label = "test$i", noise_level = Float64(i) / 10)
                    for i in 1:5
            ]
            all_scenarios = StructVector(scenarios)

            # Complete scenarios 2 and 4
            result2 = create_test_result(scenarios[2])
            result4 = create_test_result(scenarios[4])
            completed_results = StructVector([result2, result4])

            missing = find_missing_scenarios(all_scenarios, completed_results)

            @test length(missing) == 3
            @test missing[1].noise_level ≈ 0.1
            @test missing[2].noise_level ≈ 0.3
            @test missing[3].noise_level ≈ 0.5
        end
    end
end
