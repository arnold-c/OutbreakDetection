module ODStructs
export SimTimeParameters,
    EnsembleSpecification, DynamicsParameters, StateParameters

using DrWatson
@quickactivate "OutbreakDetection"

using LabelledArrays
include(srcdir("Julia/DrWatson-helpers.jl"))
include(funsdir("transmission-functions.jl"))

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

const POPULATION_N = 500_000
const LATENT_PER_DAYS = 8
const DUR_INF_DAYS = 5
const R0 = 10.0
const SIGMA = 1 / LATENT_PER_DAYS
const GAMMA = 1 / DUR_INF_DAYS
const MU = 1 / (62.5 * 365)
const BETA_MEAN = calculate_beta(R0, GAMMA, MU, 1, POPULATION_N)
const BETA_FORCE = 0.2
const EPSILON = calculate_import_rate(MU, R0, POPULATION_N)

@kwdef struct DynamicsParameters
    beta_mean::Float64 = BETA_MEAN
    beta_force::Float64 = BETA_FORCE
    sigma::Float64 = SIGMA
    gamma::Float64 = GAMMA
    mu::Float64 = MU
    epsilon::Float64 = EPSILON
    R_0::Float64 = R0
end

struct StateParameters
    init_states
    init_state_props
end

function StateParameters(;
    N = 500_00, s_prop = 0.1, e_prop = 0.01, i_prop = 0.01
)
    r_prop = 1 - (s_prop + e_prop + i_prop)
    states = @LArray [
        map(x -> Int64(round(x * N)), [s_prop, e_prop, i_prop, r_prop])..., N
    ] (:S, :E, :I, :R, :N)
    state_props = @LArray [s_prop, e_prop, i_prop, r_prop] (
        :s_prop, :i_prop, :e_prop, :r_prop
    )

    return StateParameters(
        states, state_props
    )
end

end
