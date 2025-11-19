export StateParameters

"""
    StateParameters

Parameters defining the initial state of the population.

# Fields
- `init_states::SLArray`: Initial state counts (S, E, I, R, N)
- `init_state_props::SLArray`: Initial state proportions

# Constructors
    StateParameters(; N=500_00, s_prop=0.1, e_prop=0.01, i_prop=0.01)
    StateParameters(N::Int64, init_state_props::Dict)

Creates a `StateParameters` instance with specified population size and initial
proportions in each compartment.
"""
struct StateParameters{T1 <: SLArray, T2 <: SLArray}
    init_states::T1
    init_state_props::T2
end

function StateParameters(N::Int64, init_state_props::Dict)
    return StateParameters(;
        N = N,
        s_prop = init_state_props[:s_prop],
        e_prop = init_state_props[:e_prop],
        i_prop = init_state_props[:i_prop],
    )
end

function StateParameters(;
        N = 500_00, s_prop = 0.1, e_prop = 0.01, i_prop = 0.01
    )
    r_prop = 1 - (s_prop + e_prop + i_prop)

    states = SLVector(;
        S = Int64(round(s_prop * N)),
        E = Int64(round(e_prop * N)),
        I = Int64(round(i_prop * N)),
        R = Int64(round(r_prop * N)),
        N = N,
    )
    state_props = SLVector(;
        s_prop = s_prop,
        e_prop = e_prop,
        i_prop = i_prop,
        r_prop = r_prop,
    )

    return StateParameters(
        states, state_props
    )
end
