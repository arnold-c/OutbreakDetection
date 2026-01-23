@testset "validate_all_simulations_have_outbreaks" begin
    using OutbreakDetectionCore
    using StructArrays
    using StaticArrays
    using Test
    using Try
    using Dates

    """
    Create a synthetic incidence time series with controlled outbreak patterns.

    # Arguments
    - `length::Int`: Length of the time series
    - `outbreak_periods::Vector{Tuple{Int,Int}}`: List of (start, end) indices for outbreaks
    - `baseline::Int`: Baseline incidence level
    - `outbreak_level::Int`: Incidence level during outbreaks
    """
    function create_incidence_with_outbreaks(
            length::Int,
            outbreak_periods::Vector{Tuple{Int, Int}},
            baseline::Int = 2,
            outbreak_level::Int = 10
        )
        incidence = fill(baseline, length)
        for (start_idx, end_idx) in outbreak_periods
            incidence[start_idx:end_idx] .= outbreak_level
        end
        return incidence
    end

    """
    Create a minimal SEIRRun with specified incidence pattern.
    """
    function create_seir_run(incidence::Vector{Int})
        n = length(incidence)
        # Create dummy states (not used in outbreak detection)
        states = [SVector{5, Int64}(100000, 10, 5, 0, 100015) for _ in 1:n]
        # Create dummy Reff (not used in outbreak detection)
        Reff = fill(1.0, n)
        return SEIRRun(; states = states, incidence = incidence, Reff = Reff)
    end

    """
    Create minimal ensemble and outbreak specifications for testing.

    These are just dummy specs used for error message printing - the validation
    function doesn't actually use their values, only prints them in error messages.
    """
    function create_test_specs()
        # Create a minimal ensemble spec by reading from an existing test
        # to avoid dealing with complex type constructors
        state_params = StateParameters(; N = 100_000, s_prop = 0.95, e_prop = 0.0, i_prop = 0.0)

        time_params = SimTimeParameters(; tmin = 0.0, tmax = 100.0, tstep = 1.0)

        # Use the documented API from the struct definition
        dynamics_params = TargetDiseaseDynamicsParameters(
            R_0 = 16.0,
            latent_period = Day(10),
            infectious_duration = Day(8),
            beta_force = 0.2
        )

        dynamical_noise_params = DynamicalNoiseParameters(
            R_0 = 0.1,
            latent_period = Day(1),
            infectious_duration = Day(1),
            correlation = "none",
            poisson_component = 0.0
        )

        # Create DynamicsParameterSpecification
        # Calculate derived parameters
        sigma = calculate_sigma(dynamics_params)
        gamma = calculate_gamma(dynamics_params)
        annual_births_per_k = 27.0
        population_N = state_params.init_states.N
        mu = calculate_mu(annual_births_per_k)
        epsilon = 1.0 / population_N
        beta_mean = dynamics_params.R_0 * gamma

        dynamics_spec = DynamicsParameterSpecification(
            beta_mean = beta_mean,
            beta_force = dynamics_params.beta_force,
            seasonality = dynamics_params.seasonality,
            sigma = sigma,
            gamma = gamma,
            mu = mu,
            annual_births_per_k = annual_births_per_k,
            epsilon = epsilon,
            R_0 = dynamics_params.R_0,
            min_vaccination_coverage = dynamics_params.min_vaccination_coverage,
            max_vaccination_coverage = dynamics_params.max_vaccination_coverage
        )

        ensemble_spec = EnsembleSpecification(
            "TestDisease",
            state_params,
            time_params,
            dynamics_spec,
            dynamical_noise_params,
            3
        )

        outbreak_spec = OutbreakSpecification(5, 7, 50)

        return ensemble_spec, outbreak_spec
    end

    @testset "all simulations have outbreaks" begin
        # Create 3 simulations, each with at least one outbreak
        sim1 = create_seir_run(
            create_incidence_with_outbreaks(100, [(10, 20), (50, 60)], 2, 10)
        )
        sim2 = create_seir_run(create_incidence_with_outbreaks(100, [(30, 45)], 2, 10))
        sim3 = create_seir_run(
            create_incidence_with_outbreaks(100, [(5, 15), (70, 85)], 2, 10)
        )

        ensemble_simulation = StructVector([sim1, sim2, sim3])
        ensemble_spec, outbreak_spec = create_test_specs()

        # Calculate outbreak thresholds
        outbreak_thresholds = calculate_outbreak_thresholds(
            ensemble_simulation,
            outbreak_spec
        )

        # Validate - should pass
        result = validate_all_simulations_have_outbreaks(
            outbreak_thresholds,
            ensemble_spec,
            outbreak_spec
        )

        @test Try.isok(result)
        @test Try.unwrap(result) === nothing
    end

    @testset "simulation 2 has no outbreaks" begin
        # Create 3 simulations where simulation 2 has no outbreak
        sim1 = create_seir_run(
            create_incidence_with_outbreaks(100, [(10, 20), (50, 60)], 2, 10)
        )
        # Sim2 has only baseline incidence - no outbreaks
        sim2 = create_seir_run(fill(2, 100))
        sim3 = create_seir_run(create_incidence_with_outbreaks(100, [(30, 45)], 2, 10))

        ensemble_simulation = StructVector([sim1, sim2, sim3])
        ensemble_spec, outbreak_spec = create_test_specs()

        # Calculate outbreak thresholds
        outbreak_thresholds = calculate_outbreak_thresholds(
            ensemble_simulation,
            outbreak_spec
        )

        # Validate - should fail
        result = validate_all_simulations_have_outbreaks(
            outbreak_thresholds,
            ensemble_spec,
            outbreak_spec
        )

        @test Try.iserr(result)
        error_msg = Try.unwrap_err(result)
        @test occursin("Simulation 2", error_msg)
        @test occursin("no outbreaks detected", error_msg)
        @test occursin("Total Simulations: 3", error_msg)
    end

    @testset "first simulation has no outbreaks" begin
        # Create 3 simulations where simulation 1 has no outbreak
        sim1 = create_seir_run(fill(2, 100))  # No outbreaks
        sim2 = create_seir_run(create_incidence_with_outbreaks(100, [(10, 20)], 2, 10))
        sim3 = create_seir_run(create_incidence_with_outbreaks(100, [(30, 45)], 2, 10))

        ensemble_simulation = StructVector([sim1, sim2, sim3])
        ensemble_spec, outbreak_spec = create_test_specs()

        # Calculate outbreak thresholds
        outbreak_thresholds = calculate_outbreak_thresholds(
            ensemble_simulation,
            outbreak_spec
        )

        # Validate - should fail
        result = validate_all_simulations_have_outbreaks(
            outbreak_thresholds,
            ensemble_spec,
            outbreak_spec
        )

        @test Try.iserr(result)
        error_msg = Try.unwrap_err(result)
        @test occursin("Simulation 1", error_msg)
        @test occursin("no outbreaks detected", error_msg)
    end

    @testset "last simulation has no outbreaks" begin
        # Create 3 simulations where simulation 3 has no outbreak
        sim1 = create_seir_run(create_incidence_with_outbreaks(100, [(10, 20)], 2, 10))
        sim2 = create_seir_run(create_incidence_with_outbreaks(100, [(30, 45)], 2, 10))
        sim3 = create_seir_run(fill(2, 100))  # No outbreaks

        ensemble_simulation = StructVector([sim1, sim2, sim3])
        ensemble_spec, outbreak_spec = create_test_specs()

        # Calculate outbreak thresholds
        outbreak_thresholds = calculate_outbreak_thresholds(
            ensemble_simulation,
            outbreak_spec
        )

        # Validate - should fail
        result = validate_all_simulations_have_outbreaks(
            outbreak_thresholds,
            ensemble_spec,
            outbreak_spec
        )

        @test Try.iserr(result)
        error_msg = Try.unwrap_err(result)
        @test occursin("Simulation 3", error_msg)
    end

    @testset "outbreak too short - does not qualify" begin
        # Create simulation with outbreak that's too short (only 5 days, need 7)
        sim1 = create_seir_run(create_incidence_with_outbreaks(100, [(10, 20)], 2, 10))
        sim2 = create_seir_run(
            create_incidence_with_outbreaks(100, [(30, 34)], 2, 10)
        )  # Only 5 days
        sim3 = create_seir_run(create_incidence_with_outbreaks(100, [(50, 65)], 2, 10))

        ensemble_simulation = StructVector([sim1, sim2, sim3])
        ensemble_spec, outbreak_spec = create_test_specs()

        # Calculate outbreak thresholds
        outbreak_thresholds = calculate_outbreak_thresholds(
            ensemble_simulation,
            outbreak_spec
        )

        # Validate - should fail because sim2's outbreak is too short
        result = validate_all_simulations_have_outbreaks(
            outbreak_thresholds,
            ensemble_spec,
            outbreak_spec
        )

        @test Try.iserr(result)
        error_msg = Try.unwrap_err(result)
        @test occursin("Simulation 2", error_msg)
    end

    @testset "error message includes specifications" begin
        # Create simulation with no outbreaks
        sim1 = create_seir_run(fill(2, 100))

        ensemble_simulation = StructVector([sim1])
        ensemble_spec, outbreak_spec = create_test_specs()

        # Calculate outbreak thresholds
        outbreak_thresholds = calculate_outbreak_thresholds(
            ensemble_simulation,
            outbreak_spec
        )

        # Validate - should fail
        result = validate_all_simulations_have_outbreaks(
            outbreak_thresholds,
            ensemble_spec,
            outbreak_spec
        )

        @test Try.iserr(result)
        error_msg = Try.unwrap_err(result)

        # Check that error message includes key information
        @test occursin("Ensemble Specification:", error_msg)
        @test occursin("Outbreak Specification:", error_msg)
        @test occursin("Simulation Index: 1", error_msg)
        @test occursin("Total Simulations: 1", error_msg)
        @test occursin("Please adjust either:", error_msg)
    end

    @testset "multiple simulations without outbreaks - reports first" begin
        # Create 3 simulations where 2 and 3 have no outbreaks
        sim1 = create_seir_run(create_incidence_with_outbreaks(100, [(10, 20)], 2, 10))
        sim2 = create_seir_run(fill(2, 100))  # No outbreaks
        sim3 = create_seir_run(fill(2, 100))  # No outbreaks

        ensemble_simulation = StructVector([sim1, sim2, sim3])
        ensemble_spec, outbreak_spec = create_test_specs()

        # Calculate outbreak thresholds
        outbreak_thresholds = calculate_outbreak_thresholds(
            ensemble_simulation,
            outbreak_spec
        )

        # Validate - should fail and report the first failure (sim 2)
        result = validate_all_simulations_have_outbreaks(
            outbreak_thresholds,
            ensemble_spec,
            outbreak_spec
        )

        @test Try.iserr(result)
        error_msg = Try.unwrap_err(result)
        @test occursin("Simulation 2", error_msg)
        # Should not mention simulation 3 since it fails on first error
        @test !occursin("Simulation 3", error_msg)
    end
end
