using OutbreakDetectionCore
using Test
using StaticArrays

@testset "SEIR-model.jl" verbose = true begin
    states = StateParameters(;
        N = 500_000,
        s_prop = 0.05,
        e_prop = 0.0,
        i_prop = 0.0,
    )

    target_dynamics = TargetDiseaseDynamicsParameters(
        R_0 = 16.0,
        latent_period = OutbreakDetectionCore.Dates.Day(10),
        infectious_duration = OutbreakDetectionCore.Dates.Day(8),
        beta_force = 0.2,
        min_vaccination_coverage = 0.8,
        max_vaccination_coverage = 0.8
    )
    dynamics_spec = DynamicsParameterSpecification(target_dynamics)
    dynamics_params = DynamicsParameters(dynamics_spec; seed = 1234)

    time_params = SimTimeParameters(;
        tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
    )

    seir_res = seir_mod(
        StaticArrays.SVector(states.init_states), dynamics_params, time_params; seed = 1234
    )

    @test seir_res isa SEIRRun
    @test length(seir_res.states) == time_params.tlength
    @test length(seir_res.incidence) == time_params.tlength
    @test length(seir_res.Reff) == time_params.tlength

    @test typeof(seir_res.states) <: Vector{<:StaticArrays.SVector}
    @test typeof(seir_res.incidence) <: Vector{<:Integer}
    @test typeof(seir_res.Reff) <: Vector{<:AbstractFloat}

    @test length(seir_res.states[1]) == 5
end
