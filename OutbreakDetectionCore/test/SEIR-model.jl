using OutbreakDetectionCore
using Test
using LabelledArrays: LabelledArrays
using StaticArrays

@testset "SEIR-model.jl" verbose = true begin
    states = StateParameters(;
        N = 500_000,
        s_prop = 0.05,
        e_prop = 0.00,
        i_prop = 0.00,
    )

    dynamics_params = DynamicsParameters(
        BETA_MEAN,
        BETA_FORCE,
        cos,
        SIGMA,
        GAMMA,
        MU,
        ANNUAL_BIRTHS_PER_K,
        EPSILON,
        R0,
        VACCINATION_COVERAGE,
    )

    time_params = SimTimeParameters(;
        tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
    )

    seir_res = seir_mod(
        states, dynamics_params, time_params; seed = 1234
    )

    @test length(seir_res) == 3
    seir_state_vec, seir_inc_vec, seir_beta_vec = seir_res

    @test typeof(seir_state_vec) <: Vector{<:LabelledArrays.SLArray}
    @test typeof(seir_inc_vec) <: Vector{<:StaticVector{1,<:Integer}}
    @test typeof(seir_beta_vec) <: Vector{<:AbstractFloat}

    @test length(seir_state_vec) == length(seir_inc_vec) ==
        length(
            seir_beta_vec
        ) == time_params.tlength

    @test size(seir_state_vec[1]) == 5,)
end
