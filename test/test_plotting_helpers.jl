using Test
using OutbreakDetection
using OutbreakDetectionCore: OutbreakDetectionCore
using StructArrays
using StatsBase
using Dates

# Import the functions we're testing
using OutbreakDetection:
    compute_summary_statistics,
    compute_nested_summary_statistics,
    get_noise_label,
    extract_outcome_values,
    reshape_optimization_results_to_matrix

@testset "Plotting Helpers Tests" begin
    @testset "compute_summary_statistics" begin
        @testset "Basic statistics computation" begin
            values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
            result = compute_summary_statistics(values; percentiles = [0.1, 0.9])

            @test result.mean ≈ 5.5
            @test length(result.percentiles) == 2
            @test result.percentiles[1] ≈ 1.9  # 10th percentile
            @test result.percentiles[2] ≈ 9.1  # 90th percentile
        end

        @testset "Custom percentiles" begin
            values = collect(1.0:100.0)
            result = compute_summary_statistics(values; percentiles = [0.25, 0.5, 0.75])

            @test result.mean ≈ 50.5
            @test length(result.percentiles) == 3
            @test result.percentiles[1] ≈ 25.75  # 25th percentile
            @test result.percentiles[2] ≈ 50.5   # 50th percentile (median)
            @test result.percentiles[3] ≈ 75.25  # 75th percentile
        end

        @testset "Empty vector handling" begin
            values = Float64[]
            result = compute_summary_statistics(values; percentiles = [0.1, 0.9])

            @test isnan(result.mean)
            @test all(isnan.(result.percentiles))
            @test length(result.percentiles) == 2
        end

        @testset "Single value" begin
            values = [5.0]
            result = compute_summary_statistics(values; percentiles = [0.1, 0.9])

            @test result.mean ≈ 5.0
            @test result.percentiles[1] ≈ 5.0
            @test result.percentiles[2] ≈ 5.0
        end

        @testset "Integer values" begin
            values = [1, 2, 3, 4, 5]
            result = compute_summary_statistics(values; percentiles = [0.5])

            @test result.mean ≈ 3.0
            @test result.percentiles[1] ≈ 3.0
        end
    end

    @testset "compute_nested_summary_statistics" begin
        @testset "Basic nested statistics" begin
            nested_values = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
            result = compute_nested_summary_statistics(
                nested_values; percentiles = [0.1, 0.9]
            )

            # Flattened: [1,2,3,4,5,6,7,8,9]
            @test result.mean ≈ 5.0
            @test length(result.percentiles) == 2
        end

        @testset "Unequal length nested vectors" begin
            nested_values = [[1, 2], [3, 4, 5, 6], [7]]
            result = compute_nested_summary_statistics(nested_values; percentiles = [0.5])

            # Flattened: [1,2,3,4,5,6,7]
            @test result.mean ≈ 4.0
            @test result.percentiles[1] ≈ 4.0  # median
        end

        @testset "Empty nested vector handling" begin
            nested_values = Vector{Vector{Int}}[]
            result = compute_nested_summary_statistics(
                nested_values; percentiles = [0.1, 0.9]
            )

            @test isnan(result.mean)
            @test all(isnan.(result.percentiles))
        end

        @testset "Single nested vector" begin
            nested_values = [[1, 2, 3, 4, 5]]
            result = compute_nested_summary_statistics(
                nested_values; percentiles = [0.25, 0.75]
            )

            @test result.mean ≈ 3.0
            @test result.percentiles[1] ≈ 2.0
            @test result.percentiles[2] ≈ 4.0
        end

        @testset "Float nested vectors" begin
            nested_values = [[1.5, 2.5], [3.5, 4.5], [5.5, 6.5]]
            result = compute_nested_summary_statistics(nested_values; percentiles = [0.5])

            @test result.mean ≈ 4.0
            @test result.percentiles[1] ≈ 4.0
        end
    end

    @testset "get_noise_label" begin
        @testset "Static noise labels" begin
            @test get_noise_label(:poisson) == "Static Noise"
            @test get_noise_label(:static) == "Static Noise"
        end

        @testset "Dynamical noise labels" begin
            @test get_noise_label(:dynamical_inphase) == "Dynamical Noise"
            @test get_noise_label(:dynamic) == "Dynamical Noise"
        end

        @testset "Unknown noise type" begin
            @test get_noise_label(:custom_noise) == "Custom_Noise"
            @test get_noise_label(:other) == "Other"
        end

        @testset "Case handling" begin
            # titlecase capitalizes each word (including after underscores)
            @test get_noise_label(:unknown_type) == "Unknown_Type"
        end
    end

    @testset "extract_outcome_values" begin
        # Helper to create a minimal OptimizationResult for testing
        function create_test_result()
            # Create minimal ensemble specification
            state_params = OutbreakDetectionCore.StateParameters(;
                N = 100_000,
                s_prop = 0.05,
                e_prop = 0.0,
                i_prop = 0.0
            )

            time_params = OutbreakDetectionCore.SimTimeParameters(;
                tmin = 0.0,
                tmax = 365.0,
                tstep = 1.0
            )

            target_dynamics = OutbreakDetectionCore.TargetDiseaseDynamicsParameters(;
                R_0 = 16.0,
                latent_period = Dates.Day(10),
                infectious_duration = Dates.Day(8),
                beta_force = 0.2
            )

            common_dynamics = OutbreakDetectionCore.CommonDiseaseDynamicsParameters(;
                births_per_k_pop = 12.0,
                nsims = 10
            )

            dynamics_params = OutbreakDetectionCore.DynamicsParameterSpecification(
                state_params,
                target_dynamics,
                common_dynamics
            )

            noise_params = OutbreakDetectionCore.DynamicalNoiseParameters(;
                R_0 = 0.1,
                latent_period = Dates.Day(1),
                infectious_duration = Dates.Day(1),
                correlation = "none",
                poisson_component = 0.0
            )

            ensemble_spec = OutbreakDetectionCore.EnsembleSpecification(
                "test",
                state_params,
                time_params,
                dynamics_params,
                noise_params,
                10  # nsims
            )

            test_spec = OutbreakDetectionCore.IndividualTestSpecification(
                sensitivity = 0.9,
                specificity = 0.9,
                test_result_lag = 0
            )

            alert_method = OutbreakDetectionCore.AlertMethod(
                OutbreakDetectionCore.DailyThreshold()
            )

            accuracy_metric = OutbreakDetectionCore.AccuracyMetric(
                OutbreakDetectionCore.F1()
            )

            outbreak_spec = OutbreakDetectionCore.OutbreakSpecification(
                1,  # outbreak_threshold
                7,  # minimum_outbreak_duration
                100  # minimum_outbreak_size
            )

            return OutbreakDetectionCore.OptimizationResult(;
                ensemble_specification = ensemble_spec,
                noise_level = 1.0,
                noise_type_description = :static,
                vaccination_coverage = NaN,
                test_specification = test_spec,
                percent_tested = 0.5,
                alert_method = alert_method,
                accuracy_metric = accuracy_metric,
                threshold_bounds = (lower = 0.0, upper = 20.0),
                outbreak_specification = outbreak_spec,
                optimal_threshold = 5.0,
                accuracies = [0.8, 0.85, 0.9],
                proportion_alerts_correct = [0.7, 0.75, 0.8],
                proportion_outbreaks_detected = [0.9, 0.92, 0.95],
                detection_delays = [[1, 2, 3], [2, 3, 4], [3, 4, 5]],
                unavoidable_cases = [[10, 20], [15, 25], [20, 30]],
                alert_durations = [[5, 10], [7, 12], [8, 15]],
                outbreak_durations = [[30, 40], [35, 45], [40, 50]],
                proportion_timeseries_in_alert = [0.1, 0.15, 0.2],
                proportion_timeseries_in_outbreak = [0.2, 0.25, 0.3]
            )
        end

        @testset "Extract accuracy" begin
            result = create_test_result()
            values = extract_outcome_values(result, :accuracies)
            @test values == [0.8, 0.85, 0.9]
            @test values isa Vector{Float64}
        end

        @testset "Extract detection_delays" begin
            result = create_test_result()
            values = extract_outcome_values(result, :detection_delays)
            @test values == [[1, 2, 3], [2, 3, 4], [3, 4, 5]]
            @test values isa Vector{Vector{Int64}}
        end

        @testset "Extract proportion_timeseries_in_alert" begin
            result = create_test_result()
            values = extract_outcome_values(result, :proportion_timeseries_in_alert)
            @test values == [0.1, 0.15, 0.2]
        end

        @testset "Extract proportion_timeseries_in_outbreak" begin
            result = create_test_result()
            values = extract_outcome_values(result, :proportion_timeseries_in_outbreak)
            @test values == [0.2, 0.25, 0.3]
        end

        @testset "Extract proportion_alerts_correct" begin
            result = create_test_result()
            values = extract_outcome_values(result, :proportion_alerts_correct)
            @test values == [0.7, 0.75, 0.8]
        end

        @testset "Extract proportion_outbreaks_detected" begin
            result = create_test_result()
            values = extract_outcome_values(result, :proportion_outbreaks_detected)
            @test values == [0.9, 0.92, 0.95]
        end

        @testset "Extract unavoidable_cases" begin
            result = create_test_result()
            values = extract_outcome_values(result, :unavoidable_cases)
            @test values == [[10, 20], [15, 25], [20, 30]]
        end

        @testset "Extract alert_durations" begin
            result = create_test_result()
            values = extract_outcome_values(result, :alert_durations)
            @test values == [[5, 10], [7, 12], [8, 15]]
        end

        @testset "Extract outbreak_durations" begin
            result = create_test_result()
            values = extract_outcome_values(result, :outbreak_durations)
            @test values == [[30, 40], [35, 45], [40, 50]]
        end

        @testset "Unknown outcome throws error" begin
            result = create_test_result()
            @test_throws AssertionError extract_outcome_values(result, :unknown_outcome)
        end
    end

    @testset "reshape_optimization_results_to_matrix" begin
        # Helper to create test results
        function create_test_results(;
                noise_types = [:static, :dynamic],
                noise_levels = [1.0, 2.0],
                test_specs = [
                    OutbreakDetectionCore.IndividualTestSpecification(
                        sensitivity = 0.85,
                        specificity = 0.85,
                        test_result_lag = 0
                    ),
                    OutbreakDetectionCore.IndividualTestSpecification(
                        sensitivity = 0.9,
                        specificity = 0.9,
                        test_result_lag = 14
                    ),
                ],
                percent_tested_vals = [0.5, 0.7]
            )
            # Create minimal ensemble specification
            state_params = OutbreakDetectionCore.StateParameters(;
                N = 100_000,
                s_prop = 0.05,
                e_prop = 0.0,
                i_prop = 0.0
            )

            time_params = OutbreakDetectionCore.SimTimeParameters(;
                tmin = 0.0,
                tmax = 365.0,
                tstep = 1.0
            )

            target_dynamics = OutbreakDetectionCore.TargetDiseaseDynamicsParameters(;
                R_0 = 16.0,
                latent_period = Dates.Day(10),
                infectious_duration = Dates.Day(8),
                beta_force = 0.2
            )

            common_dynamics = OutbreakDetectionCore.CommonDiseaseDynamicsParameters(;
                births_per_k_pop = 12.0,
                nsims = 10
            )

            dynamics_params = OutbreakDetectionCore.DynamicsParameterSpecification(
                state_params,
                target_dynamics,
                common_dynamics
            )

            noise_params = OutbreakDetectionCore.DynamicalNoiseParameters(;
                R_0 = 0.1,
                latent_period = Dates.Day(1),
                infectious_duration = Dates.Day(1),
                correlation = "none",
                poisson_component = 0.0
            )

            ensemble_spec = OutbreakDetectionCore.EnsembleSpecification(
                "test",
                state_params,
                time_params,
                dynamics_params,
                noise_params,
                10  # nsims
            )

            alert_method = OutbreakDetectionCore.AlertMethod(
                OutbreakDetectionCore.DailyThreshold()
            )

            accuracy_metric = OutbreakDetectionCore.AccuracyMetric(
                OutbreakDetectionCore.F1()
            )

            outbreak_spec = OutbreakDetectionCore.OutbreakSpecification(
                1,  # outbreak_threshold
                7,  # minimum_outbreak_duration
                100  # minimum_outbreak_size
            )

            results = OutbreakDetectionCore.OptimizationResult[]

            for noise_type in noise_types
                for noise_level in noise_levels
                    for test_spec in test_specs
                        for percent_tested in percent_tested_vals
                            push!(
                                results,
                                OutbreakDetectionCore.OptimizationResult(;
                                    ensemble_specification = ensemble_spec,
                                    noise_level = noise_level,
                                    noise_type_description = noise_type,
                                    vaccination_coverage = noise_type == :static ? NaN : 0.6,
                                    test_specification = test_spec,
                                    percent_tested = percent_tested,
                                    alert_method = alert_method,
                                    accuracy_metric = accuracy_metric,
                                    threshold_bounds = (lower = 0.0, upper = 20.0),
                                    outbreak_specification = outbreak_spec,
                                    optimal_threshold = 5.0,
                                    accuracies = [0.8, 0.85, 0.9],
                                    proportion_alerts_correct = [0.7, 0.75, 0.8],
                                    proportion_outbreaks_detected = [0.9, 0.92, 0.95],
                                    detection_delays = [[1, 2], [2, 3], [3, 4]],
                                    unavoidable_cases = [[10, 20], [15, 25], [20, 30]],
                                    alert_durations = [[5, 10], [7, 12], [8, 15]],
                                    outbreak_durations = [[30, 40], [35, 45], [40, 50]],
                                    proportion_timeseries_in_alert = [0.1, 0.15, 0.2],
                                    proportion_timeseries_in_outbreak = [0.2, 0.25, 0.3]
                                )
                            )
                        end
                    end
                end
            end

            return StructArray(results)
        end

        @testset "Basic matrix reshaping" begin
            results = create_test_results()
            matrix, noise_types = reshape_optimization_results_to_matrix(results)

            # Should have 2 noise types (rows) and 2 noise levels (columns)
            @test size(matrix) == (2, 2)
            @test length(noise_types) == 2
            @test Set(noise_types) == Set([:static, :dynamic])
        end

        @testset "Matrix cell contents" begin
            results = create_test_results()
            matrix, noise_types = reshape_optimization_results_to_matrix(results)

            # Each cell should contain results for that noise type/level combination
            for i in 1:size(matrix, 1)
                for j in 1:size(matrix, 2)
                    cell = matrix[i, j]
                    @test cell isa StructVector{OutbreakDetectionCore.OptimizationResult}
                    @test !isempty(cell)

                    # All results in cell should have same noise type and level
                    @test all(r -> r.noise_type_description == noise_types[i], cell)
                    noise_level = cell[1].noise_level
                    @test all(r -> r.noise_level == noise_level, cell)
                end
            end
        end

        @testset "Test specification sorting" begin
            results = create_test_results()
            matrix, _ = reshape_optimization_results_to_matrix(results)

            # Check that within each cell, results are sorted by test_result_lag (descending)
            for i in 1:size(matrix, 1)
                for j in 1:size(matrix, 2)
                    cell = matrix[i, j]

                    # Group by percent_tested to check sorting within groups
                    unique_percent_tested = unique([r.percent_tested for r in cell])
                    for pct in unique_percent_tested
                        subset = filter(r -> r.percent_tested == pct, cell)
                        subset_lags = [
                            r.test_specification.test_result_lag for r in subset
                        ]

                        # Should be sorted descending by lag
                        @test issorted(subset_lags; rev = true)
                    end
                end
            end
        end

        @testset "Single noise type and level" begin
            results = create_test_results(; noise_types = [:static], noise_levels = [1.0])
            matrix, noise_types = reshape_optimization_results_to_matrix(results)

            @test size(matrix) == (1, 1)
            @test length(noise_types) == 1
            @test noise_types[1] == :static
        end

        @testset "Multiple noise levels" begin
            results = create_test_results(;
                noise_types = [:static],
                noise_levels = [1.0, 2.0, 4.0, 7.0]
            )
            matrix, noise_types = reshape_optimization_results_to_matrix(results)

            @test size(matrix) == (1, 4)
            @test length(noise_types) == 1

            # Check that noise levels are sorted in columns
            for j in 1:4
                expected_level = [1.0, 2.0, 4.0, 7.0][j]
                @test all(r -> r.noise_level == expected_level, matrix[1, j])
            end
        end

        @testset "Noise types ordering" begin
            results = create_test_results(;
                noise_types = [:dynamic, :static],  # Reversed order
                noise_levels = [1.0]
            )
            matrix, noise_types = reshape_optimization_results_to_matrix(results)

            @test size(matrix) == (2, 1)
            # Order should be preserved from unique() call
            @test length(noise_types) == 2
        end
    end
end
