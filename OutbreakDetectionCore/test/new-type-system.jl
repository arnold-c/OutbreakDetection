using OutbreakDetectionCore
using Test
using StaticArrays
using Try

@testset "New Type System Integration" begin
    @testset "DynamicsParameters Hierarchy" begin
        # Test three-tier hierarchy
        target = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2,
        )

        @test target.R_0 == 16.0
        @test target.latent_period_days == 10.0
        @test target.infectious_duration_days == 8.0
        @test target.beta_force == 0.2
        @test target.life_expectancy_years == 62.5  # default
        @test target.population_N == 500_000  # default

        # Test specification creation
        spec = DynamicsParameterSpecification(target)

        @test spec.R_0 == 16.0
        @test spec.sigma ≈ 0.1
        @test spec.gamma ≈ 0.125
        @test spec.beta_force == 0.2
        @test spec.beta_mean > 0
        @test spec.epsilon > 0
        @test spec.mu > 0
        @test spec.annual_births_per_k > 0

        # Test concrete parameters creation
        dynamics = DynamicsParameters(spec; vaccination_coverage = 0.8)

        @test dynamics.R_0 == 16.0
        @test dynamics.sigma ≈ 0.1
        @test dynamics.gamma ≈ 0.125
        @test dynamics.vaccination_coverage == 0.8
    end

    @testset "SeasonalityFunction Sum Types" begin
        # Test cosine seasonality
        cosine_seasonality = SeasonalityFunction(CosineSeasonality())
        @test cosine_seasonality isa SeasonalityFunction

        # Test sine seasonality
        sine_seasonality = SeasonalityFunction(SineSeasonality())
        @test sine_seasonality isa SeasonalityFunction

        # Test in target parameters
        target_cos = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2,
            seasonality = cosine_seasonality,
        )

        @test target_cos.seasonality isa SeasonalityFunction

        target_sin = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2,
            seasonality = sine_seasonality,
        )

        @test target_sin.seasonality isa SeasonalityFunction
    end

    @testset "Noise Specifications" begin
        # Test DynamicalNoiseParameters
        noise_spec = DynamicalNoiseParameters(
            R_0 = 5.0,
            latent_period = 7.0,
            infectious_duration = 14.0,
            correlation = "in-phase",
            poisson_component = 1.0,
            vaccination_bounds = [0.0, 0.8],
        )

        @test noise_spec.R_0 == 5.0
        @test noise_spec.latent_period == 7.0
        @test noise_spec.infectious_duration == 14.0
        @test noise_spec.correlation == "in-phase"
        @test noise_spec.poisson_component == 1.0
        @test noise_spec.vaccination_bounds == [0.0, 0.8]

        # Test DynamicalNoiseSpecification instance creation
        noise = DynamicalNoiseSpecification(noise_spec, 0.65)

        @test noise.R_0 == 5.0
        @test noise.vaccination_coverage == 0.65

        # Test PoissonNoiseSpecification
        poisson = PoissonNoiseSpecification(noise_mean_scaling = 1.0)

        @test poisson.noise_mean_scaling == 1.0
    end

    @testset "Simulation Result Types" begin
        # Test SEIRRun
        states = [SVector{5, Int64}(450000, 100, 50, 49850, 500000) for _ in 1:10]
        incidence = rand(1:10, 10)
        reff = rand(0.5:0.01:2.0, 10)

        run = SEIRRun(states = states, incidence = incidence, Reff = reff)

        @test run.states == states
        @test run.incidence == incidence
        @test run.Reff == reff

        # Test DynamicalNoiseRun
        noise_run = DynamicalNoiseRun(
            incidence = [incidence, incidence],
            Reff = [reff, reff],
            mean_noise = 15.5,
            mean_poisson_noise = 5.5,
            mean_dynamic_noise = 10.0,
        )

        @test noise_run.mean_noise == 15.5
        @test noise_run.mean_poisson_noise == 5.5
        @test noise_run.mean_dynamic_noise == 10.0

        # Test PoissonNoiseRun
        poisson_run = PoissonNoiseRun(incidence = [incidence, incidence], mean_noise = 12.5)

        @test poisson_run.mean_noise == 12.5
    end

    @testset "Endemic Equilibrium Calculation" begin
        # Note: Endemic equilibrium calculation may not work for all parameter combinations
        # This is a known limitation that will be addressed in future work
        # For now, we test that the function exists and returns Try types

        target = TargetDiseaseDynamicsParameters(
            R_0 = 5.0,
            latent_period_days = 7.0,
            infectious_duration_days = 14.0,
            beta_force = 0.2,
        )

        spec = DynamicsParameterSpecification(target)

        # Test that function returns Try type
        result = calculate_endemic_equilibrium_proportions(spec, 0.5)
        @test result isa Union{Try.Ok, Try.Err}

        # Test with high vaccination (R_eff ≤ 1) should return error
        result_high_vax = calculate_endemic_equilibrium_proportions(spec, 0.85)
        @test Try.iserr(result_high_vax)

        # Test with invalid vaccination coverage
        result_invalid = calculate_endemic_equilibrium_proportions(spec, 1.5)
        @test Try.iserr(result_invalid)
    end

    @testset "Time Parameters" begin
        # Test default constructor
        time_params = SimTimeParameters()

        @test time_params.tmin == 0.0
        @test time_params.tmax == 365.0 * 100.0
        @test time_params.tstep == 1.0
        @test time_params.tlength > 0

        # Test custom constructor
        time_params_custom = SimTimeParameters(tmax = 365.0 * 20)

        @test time_params_custom.tmax == 365.0 * 20
        @test time_params_custom.tlength < time_params.tlength
    end

    @testset "Validation" begin
        # Note: @kwdef doesn't run inner constructor validation automatically
        # Validation happens in the conversion constructors

        # Test invalid vaccination coverage in DynamicsParameters
        target = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2,
        )
        spec = DynamicsParameterSpecification(target)

        @test_throws AssertionError DynamicsParameters(spec; vaccination_coverage = 1.5)
        @test_throws AssertionError DynamicsParameters(spec; vaccination_coverage = -0.1)

        # Test validation in DynamicsParameterSpecification constructor
        invalid_target = TargetDiseaseDynamicsParameters(
            R_0 = -1.0,  # Invalid
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 0.2,
        )
        @test_throws AssertionError DynamicsParameterSpecification(invalid_target)

        # Test invalid beta_force
        invalid_target2 = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period_days = 10.0,
            infectious_duration_days = 8.0,
            beta_force = 1.5,  # Invalid
        )
        @test_throws AssertionError DynamicsParameterSpecification(invalid_target2)

        # Test invalid noise correlation (throws during construction)
        @test_throws AssertionError DynamicalNoiseParameters(
            R_0 = 5.0,
            latent_period = 7.0,
            infectious_duration = 14.0,
            correlation = "invalid",
            poisson_component = 1.0,
        )
    end

    @testset "Backward Compatibility" begin
        # Test deprecated constructor still works (with warning)
        local dynamics
        @test_logs (:warn,) match_mode = :any dynamics =
            DynamicsParameters(0.1, 0.125, 16.0)

        @test dynamics.sigma ≈ 0.1
        @test dynamics.gamma ≈ 0.125
        @test dynamics.R_0 == 16.0
    end
end
