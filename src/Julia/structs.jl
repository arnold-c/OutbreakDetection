using LabelledArrays

struct EnsembleTransmissionParameters
    R_0::Float64
    sigma::Float64
    gamma::Float64
end

struct SimTimeParameters
    tmin
    tmax
    tstep
    trange
    tlength
end

function SimTimeParameters(; tmin = 0.0, tmax = 365.0 * 100.0, tstep = 1.0)
    return SimTimeParameters(
        tmin, tmax, tstep, tmin:tstep:tmax, length(tmin:tstep:tmax)
    )
end

struct EnsembleSpecification
    modeltypes::Tuple
    N::Int64
    Rinit_prop::Float64
    nsims::Int64
    births_per_k::Int64
    beta_force::Float64
    time_parameters::SimTimeParameters
end

@kwdef struct DynamicsParameters1
    beta_mean::Float64 = 4.00e-8
    beta_force::Float64 = 0.2
    sigma::Float64 = 0.125
    gamma::Float64 = 0.2
    mu::Float64 = 4.38e-5
    epsilon::Float64 = 5.91e-7
    R_0::Float64 = 10.0
end

DynamicsParameters = DynamicsParameters1

struct StateParameters2
    init_states
    init_state_props
end

function StateParameters(;
    N = 500_00, s_prop = 0.1, e_prop = 0.01, i_prop = 0.01
)
    r_prop = 1 - (s_prop + e_prop + i_prop)
    states = @LArray [map(x -> Int64(round(x * N)), [s_prop, e_prop, i_prop, r_prop])..., N] (:S, :E, :I, :R, :N)
    state_props = @LArray [s_prop, e_prop, i_prop, r_prop] (:s_prop, :i_prop, :e_prop, :r_prop)

    return StateParameters2(
        states, state_props
    )
end
