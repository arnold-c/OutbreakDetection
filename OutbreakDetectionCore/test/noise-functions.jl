@testset "noise-functions.jl" begin
    using OutbreakDetectionCore
    using StaticArrays
    using NaNMath: NaNMath
    using Statistics
    using StructArrays

    # New API: DynamicalNoiseSpecification
    dynamical_noise_spec = DynamicalNoiseSpecification(
        R_0 = 5.0,
        latent_period = 7.0,
        duration_infection = 14.0,
        correlation = "in-phase",
        poisson_component = 0.15,
        vaccination_bounds = [0.0, 0.8]
    )
    # Create concrete instance for simulation
    dynamical_noise = DynamicalNoise(dynamical_noise_spec, 0.0) # 0.0 vaccination for test

    state_parameters = StateParameters(
        500_000,
        Dict(
            :s_prop => 0.05,
            :e_prop => 0.0,
            :i_prop => 0.001,
            :r_prop => 0.949,
        ),
    )

    # New API: DynamicsParameters
    target_dynamics = TargetDiseaseDynamicsParameters(
        R_0 = 16.0,
        latent_period_days = 10.0,
        infectious_duration_days = 8.0,
        beta_force = 0.2,
        population_N = 500_000
    )
    dynamics_spec = DynamicsParameterSpecification(target_dynamics)
    dynamics_parameters = DynamicsParameters(dynamics_spec; vaccination_coverage = 0.8)

    time_parameters = SimTimeParameters(;
        tmin = 0.0,
        tmax = 365.0 * 100.0,
        tstep = 1.0,
    )
    nsims = 100
    seed = 1234
    ensemble_spec = EnsembleSpecification(
        ("seasonal-infectivity-import", "tau-leaping"),
        state_parameters,
        dynamics_parameters,
        time_parameters,
        nsims,
    )

    # Run SEIR simulation using new API
    # We need SEIR results for Poisson noise generation
    seir_results = StructVector{SEIRRun}(undef, nsims)

    # Helper to run single simulation (since we don't have ensemble run function exposed here yet)
    for sim in 1:nsims
        run_seed = seed + (sim - 1)
        seir_results[sim] = seir_mod(
            StaticArrays.SVector(state_parameters.init_states),
            dynamics_parameters,
            time_parameters;
            seed = run_seed
        )
    end

    # Test Dynamical Noise
    # Create noise dynamics parameters for the test
    noise_dynamics_params = create_noise_dynamics_parameters(
        dynamical_noise,
        dynamics_spec,
        state_parameters.init_states.N
    )

    noise_run = create_noise_vecs(
        dynamical_noise,
        ensemble_spec,
        noise_dynamics_params;
        seed = seed
    )

    @test noise_run isa NoiseRun
    @test noise_run.mean_noise ≈ noise_run.mean_dynamic_noise + noise_run.mean_poisson_noise

    # Test Poisson Noise
    poisson_noise = PoissonNoise(0.15) # 0.15 scaling

    poisson_run = create_noise_vecs(
        poisson_noise,
        ensemble_spec,
        seir_results;
        seed = seed
    )

    @test poisson_run isa NoiseRun
    @test poisson_run.mean_noise > 0.0
end
