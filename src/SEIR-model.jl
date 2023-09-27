# module SEIRModel
#
# export calculate_beta_amp, seir_mod, seir_mod!, seir_mod_loop!
#
"""
This is a simulation of an SIR model that uses Tau-leaping, with commuter
imports. All jumps are manually defined.
"""

using DifferentialEquations
using Statistics
using Distributions
using Random
using UnPack
using LoopVectorization
using StaticArrays

function seir_static_mod(
    states,
    dynamics_params,
    time_params;
    type = "stoch",
    seed = 1234
)
    Random.seed!(seed)

    state_vec = Vector{typeof(states)}(undef, time_params.tlength)
    change_vec = similar(state_vec)
    jump_vec = Vector{SVector{10,Float64}}(undef, time_params.tlength)

    beta_vec = map(
        t -> calculate_beta_amp(
            dynamics_params.beta_mean, dynamics_params.beta_force, t
        ),
        time_params.trange,
    )

    state_vec[1] = states
    change_vec[1] = 0.0
    jump_vec[1] = 0.0
    for i in 1:(length(time_params.trange) - 1)
        state_vec[i + 1], change_vec[i + 1], jump_vec[i + 1] = seir_static_mod_loop(
            beta_vec[i],
            states[i],
            dynamics_params,
            time_params;
            type = type
        )
    end
    return state_vec, change_vec, jump_vec, beta_vec
end

function seir_static_mod_loop(
    beta_t, states, dynamics_params, time_params; type = "stoch"
)
    S, E, I, R, N = states

    # Calculate the rates of each event
    infec_rate = beta_t * S * I                   # Contact: S -> E
    latent_rate = dynamics_params.sigma * E        # E -> I
    recov_rate = dynamics_params.gamma * I        # I -> R
    susc_birth_rate =
        dynamics_params.mu * (1 - dynamics_params.vaccination_coverage) * N           # Birth -> S
    S_death_rate = dynamics_params.mu * S           # S -> death
    E_death_rate = dynamics_params.mu * E           # E -> death
    I_death_rate = dynamics_params.mu * I           # I -> death
    R_death_rate = dynamics_params.mu * R           # R -> death
    import_rate = dynamics_params.epsilon * N / dynamics_params.R_0        # Import: S -> E
    vacc_birth_rate =
        dynamics_params.mu * dynamics_params.vaccination_coverage * N    # Birth -> R

    rates = @SVector [
        infec_rate,
        latent_rate,
        recov_rate,
        susc_birth_rate,
        S_death_rate,
        E_death_rate,
        I_death_rate,
        R_death_rate,
        import_rate,
        vacc_birth_rate,
    ]

    # Calculate the number of jumps for each event
    if type == "stoch"
        jumps = map(r -> rand(Poisson(r * time_params.tstep)), rates)
    elseif type == "det"
        jumps = map(r -> r * time_params.tstep, rates)
    else
        return ("Type must be stoch or det")
    end

    # Calculate the change in each state
    S_change = jumps[4] - jumps[1] - jumps[5] - jumps[9]
    E_change = jumps[1] - jumps[2] - jumps[6] + jumps[9]
    I_change = jumps[2] - jumps[3] - jumps[7]
    R_change = jumps[3] - jumps[8] + jumps[10]

    # Check that the change in each state does not result in a negative state,
    # and if it is, set the change to the negative of the current state (i.e.
    # set the new state value to 0)
    for (state, change) in
        zip(states, @SVector [S_change, E_change, I_change, R_change])
        if state + change < 0
            change = -state
        end
    end

    N_change = sum(@SVector [S_change, E_change, I_change, R_change])

    changes = @SVector [S_change, E_change, I_change, R_change]

    new_S = S + S_change
    new_E = E + E_change
    new_I = I + I_change
    new_R = R + R_change
    new_N = N + N_change

    new_states = @SVector [new_S, new_E, new_I, new_R, new_N]

    return new_states, changes, jumps
end

"""
    seir_mod(states, dynamics_params, trange; tstep, type = "stoch")

The in-place function to run the SEIR model with a vaccinations going directly to the R compartment and produce the transmission rate array.
"""
function seir_mod(
    states,
    dynamics_params,
    time_params;
    type = "stoch",
    seed = 1234
)
    Random.seed!(seed)

    state_arr = zeros(Float64, time_params.tlength, size(states, 1))
    change_arr = similar(state_arr)
    jump_arr = zeros(
        Float64, time_params.tlength, (size(states, 1) - 1) * 2 + 2
    )
    beta_arr = zeros(Float64, time_params.tlength)

    seir_mod!(
        state_arr,
        change_arr,
        jump_arr,
        beta_arr,
        states,
        Vector{Float64}(undef, size(jump_arr, 2)),
        dynamics_params,
        time_params;
        type = type,
        seed = seed,
    )
    return state_arr, change_arr, jump_arr, beta_arr
end

"""
    seir_mod!(state_arr, change_arr, jump_arr, beta_arr, states, dynamics_params, trange; tstep, type = "stoch")

The in-place function to run the SEIR model and produce the transmission rate array.
"""
function seir_mod!(
    state_arr,
    change_arr,
    jump_arr,
    beta_arr,
    states,
    rates,
    dynamics_params,
    time_params;
    type = "stoch",
    seed = 1234,
)
    Random.seed!(seed)

    @. beta_arr = calculate_beta_amp(
        dynamics_params.beta_mean,
        dynamics_params.beta_force,
        time_params.trange,
    )

    for i in eachindex(time_params.trange)
        if i == 1
            state_arr[i, :] .= states
            continue
        end

        seir_mod_loop!(
            state_arr,
            change_arr,
            jump_arr,
            beta_arr,
            rates,
            i,
            dynamics_params,
            time_params;
            type = type,
        )
    end

    return nothing
end

"""
    seir_mod_loop!(state_arr, change_arr, jump_arr, j, params, t, tstep; type = type)

The inner loop that is called by `seir_mod!()` function.
"""
function seir_mod_loop!(
    state_arr,
    change_arr,
    jump_arr,
    beta_arr,
    rates,
    i,
    dynamics_params,
    time_params;
    type = type,
)
    # TODO: Benchmak StaticArrays implementation as potentially much faster.
    # Would need to use permutedims(reshape(reinterperate(Float64, SVector), (...), (...))
    # to get it into an array that could be used later on.
    # Create views of the state variables for easier use
    @views S = state_arr[i - 1, 1]
    @views E = state_arr[i - 1, 2]
    @views I = state_arr[i - 1, 3]
    @views R = state_arr[i - 1, 4]
    @views N = state_arr[i - 1, 5]
    @views beta_t = beta_arr[i]

    # Calculate the rates of each event
    rates[1] = beta_t * S * I                   # Contact: S -> E
    rates[2] = dynamics_params.sigma * E        # E -> I
    rates[3] = dynamics_params.gamma * I        # I -> R
    rates[4] =
        dynamics_params.mu * (1 - dynamics_params.vaccination_coverage) * N           # Birth -> S
    rates[5] = dynamics_params.mu * S           # S -> death
    rates[6] = dynamics_params.mu * E           # E -> death
    rates[7] = dynamics_params.mu * I           # I -> death
    rates[8] = dynamics_params.mu * R           # R -> death
    rates[9] = dynamics_params.epsilon * N / dynamics_params.R_0        # Import: S -> E
    rates[10] = dynamics_params.mu * dynamics_params.vaccination_coverage * N    # Birth -> R

    # Calculate the number of jumps for each event
    if type == "stoch"
        @simd for r in eachindex(rates)
            jump_arr[i, r] = rand(Poisson(rates[r] * time_params.tstep))
        end
    elseif type == "det"
        @simd for r in eachindex(rates)
            jump_arr[i, r] = rates[r] * time_params.tstep
        end
    else
        return ("Type must be stoch or det")
    end

    # Calculate the change in each state
    @views change_arr[i, 1] =
        jump_arr[i, 4] - jump_arr[i, 1] - jump_arr[i, 5] - jump_arr[i, 9]
    @views change_arr[i, 2] =
        jump_arr[i, 1] - jump_arr[i, 2] - jump_arr[i, 6] + jump_arr[i, 9]
    @views change_arr[i, 3] = jump_arr[i, 2] - jump_arr[i, 3] - jump_arr[i, 7]
    @views change_arr[i, 4] = jump_arr[i, 3] - jump_arr[i, 8] + jump_arr[i, 10]

    # Check that the change in each state does not result in a negative state,
    # and if it is, set the change to the negative of the current state (i.e.
    # set the new state value to 0)
    for state in axes(state_arr, 2)
        @views if state_arr[i - 1, state] + change_arr[i, state] < 0
            @views change_arr[i, state] = -state_arr[i - 1, state]
        end
    end

    @views change_arr[i, 5] = sum(change_arr[i, 1:4])

    @views previous_states = state_arr[i - 1, :]
    @views current_changes = change_arr[i, :]
    @. state_arr[i, :] = previous_states + current_changes

    return nothing
end

# end
