@testset "noise-functions.jl" begin
    using OutbreakDetectionCore
    using StaticArrays
    using NaNMath: NaNMath

    dynamical_noise_spec = DynamicalNoiseSpecification(
        "dynamical",
        5.0,
        7,
        14,
        "in-phase",
        0.15,
        0.0,
    )

    state_parameters = StateParameters(
        500_000,
        Dict(
            :s_prop => 0.05,
            :e_prop => 0.0,
            :i_prop => 0.0,
            :r_prop => 0.95,
        ),
    )
    dynamics_parameters = DynamicsParameters(
        500_000,
        27,
        0.2,
        SIGMA,
        GAMMA,
        16.0,
        0.8,
    )
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

    ensemble_seir_vecs = Array{typeof(state_parameters.init_states), 2}(
        undef,
        time_parameters.tlength,
        nsims,
    )

    ensemble_inc_vecs = Array{typeof(StaticArrays.SVector(0)), 2}(
        undef,
        time_parameters.tlength,
        nsims,
    )

    ensemble_beta_arr = zeros(Float64, time_parameters.tlength)

    for sim in axes(ensemble_inc_vecs, 2)
        run_seed = seed + (sim - 1)

        seir_mod!(
            @view(ensemble_seir_vecs[:, sim]),
            @view(ensemble_inc_vecs[:, sim]),
            ensemble_beta_arr,
            state_parameters.init_states,
            dynamics_parameters,
            time_parameters;
            seed = run_seed,
        )
    end

    ensemble_inc_arr = zeros(
        Int64, size(ensemble_inc_vecs, 1), size(ensemble_inc_vecs, 2)
    )

    for sim in axes(ensemble_inc_vecs, 2)
        convert_svec_to_matrix!(
            @view(ensemble_inc_arr[:, sim]),
            @view(ensemble_inc_vecs[:, sim])
        )
    end

    dynamical_noise_array, dynamical_noise_means = create_noise_arr(
        dynamical_noise_spec,
        ensemble_inc_arr;
        ensemble_specification = ensemble_spec,
        seed = seed,
    )

    @test isapprox(
        dynamical_noise_means.mean_poisson_noise +
            dynamical_noise_means.mean_rubella_noise,
        mean(dynamical_noise_array),
        atol = 1.0e-2,
    )

    @test isequal(
        dynamical_noise_means.mean_poisson_noise +
            dynamical_noise_means.mean_rubella_noise,
        dynamical_noise_means.mean_noise,
    )

    @test isapprox(
        dynamical_noise_means.poisson_noise_prop,
        dynamical_noise_spec.noise_mean_scaling,
        atol = 1.0e-2,
    )

    poisson_noise_spec = PoissonNoiseSpecification(
        "poisson",
        0.15,
    )

    poisson_noise_array, poisson_noise_means = create_noise_arr(
        poisson_noise_spec,
        ensemble_inc_arr;
        seed = seed,
    )

    @test isapprox(
        poisson_noise_means.mean_poisson_noise,
        mean(ensemble_inc_arr * poisson_noise_spec.noise_mean_scaling),
        atol = 1.0e-2,
    )

    @test isequal(poisson_noise_means.mean_noise, mean(poisson_noise_array))
end
