using OutbreakDetectionCore
using Test
using StructArrays
using AutoHashEquals

@testset "Optimization Wrapper" begin
    @testset "ScenarioParameters" begin
        target = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2,
        )

        noise = DynamicalNoiseSpecification(
            R_0 = 5.0,
            latent_period = 7.0,
            duration_infection = 14.0,
            correlation = "in-phase",
            poisson_component = 1.0,
            vaccination_bounds = [0.0, 0.8],
        )

        test_spec = IndividualTestSpecification(1.0, 0.8, 0)

        scenario1 = ScenarioParameters(target, noise, test_spec, 0.5, 0.01)

        scenario2 = ScenarioParameters(target, noise, test_spec, 0.5, 0.01)

        # Test AutoHashEquals
        @test scenario1 == scenario2
        @test hash(scenario1) == hash(scenario2)

        # Test inequality
        scenario3 = ScenarioParameters(target, noise, test_spec, 0.7, 0.01)  # Different percent_clinic_tested

        @test scenario1 != scenario3
        @test hash(scenario1) != hash(scenario3)
    end

    @testset "OptimizationResult StructVector" begin
        # Create a scenario for testing
        target = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2,
        )

        noise = DynamicalNoiseSpecification(
            R_0 = 5.0,
            latent_period = 7.0,
            duration_infection = 14.0,
            correlation = "in-phase",
            poisson_component = 1.0,
            vaccination_bounds = [0.0, 0.8],
        )

        test_spec = IndividualTestSpecification(1.0, 0.8, 0)

        scenario = ScenarioParameters(target, noise, test_spec, 0.5, 0.01)

        # Create mock results
        results = [
            OptimizationResult(scenario, 0.05, 0.95, 0.92, 0.98, 14.5, 0.88) for i in 1:10
        ]

        # Convert to StructVector
        results_sv = StructVector(results)

        @test results_sv isa StructVector
        @test length(results_sv) == 10

        # Test column access
        @test results_sv.accuracy isa Vector{Float64}
        @test all(results_sv.accuracy .== 0.95)

        # Test filtering
        high_acc = filter(r -> r.accuracy > 0.9, results_sv)
        @test length(high_acc) == 10
    end

    @testset "Scenario Grid Creation" begin
        base_target = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2,
        )

        base_noise = DynamicalNoiseSpecification(
            R_0 = 5.0,
            latent_period = 7.0,
            duration_infection = 14.0,
            correlation = "in-phase",
            poisson_component = 1.0,
            vaccination_bounds = [0.0, 0.8],
        )

        base_test = IndividualTestSpecification(1.0, 0.8, 0)

        scenarios = create_scenario_grid(
            R_0_values = [12.0, 16.0],
            noise_levels = ["low", "high"],
            test_sensitivities = [0.9, 1.0],
            percent_tested_values = [0.5],
            outbreak_thresholds = [0.01, 0.02],
            base_target = base_target,
            base_noise = base_noise,
            base_test = base_test,
        )

        # Should have 2 × 2 × 2 × 1 × 2 = 16 scenarios
        @test length(scenarios) == 16

        # Test uniqueness (AutoHashEquals)
        @test length(unique(scenarios)) == 16

        # Test that R_0 values are correct
        R0_values = [s.target_dynamics.R_0 for s in scenarios]
        @test 12.0 in R0_values
        @test 16.0 in R0_values

        # Test that test sensitivities are correct
        sensitivities = [s.test_spec.sensitivity for s in scenarios]
        @test 0.9 in sensitivities
        @test 1.0 in sensitivities
    end

    @testset "Noise Level Mapping" begin
        base_noise = DynamicalNoiseSpecification(
            R_0 = 5.0,
            latent_period = 7.0,
            duration_infection = 14.0,
            correlation = "in-phase",
            poisson_component = 1.0,
            vaccination_bounds = [0.0, 0.8],
        )

        # Test low noise level
        low_noise = OutbreakDetectionCore._create_noise_for_level("low", base_noise)
        @test low_noise.vaccination_bounds == [0.0, 0.5]

        # Test medium noise level
        medium_noise = OutbreakDetectionCore._create_noise_for_level("medium", base_noise)
        @test medium_noise.vaccination_bounds == [0.3, 0.7]

        # Test high noise level
        high_noise = OutbreakDetectionCore._create_noise_for_level("high", base_noise)
        @test high_noise.vaccination_bounds == [0.5, 0.9]

        # Test error on unknown level
        @test_throws ErrorException OutbreakDetectionCore._create_noise_for_level(
            "unknown", base_noise
        )
    end
end
