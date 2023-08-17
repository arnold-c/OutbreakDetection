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

function SimTimeParameters(; tmin, tmax, tstep)
    return SimTimeParameters(tmin, tmax, tstep, tmin:tstep:tmax, length(tmin:tstep:tmax))
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
    beta_mean::Float64
    beta_force::Float64 = 0.2
    sigma::Float64 = 0.125
    gamma::Float64 = 0.2
    mu::Float64 = 4.38e-5
    epsilon::Float64 = 5.91e-7
    R_0::Float64 = 10.0
end

DynamicsParameters = DynamicsParameters1
