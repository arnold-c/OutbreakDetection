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


function seir_wide_mod(
    states,
    dynamics_params,
    time_params;
    retbetaarr = false,
    type = "stoch",
    seed = 1234,
)
    Random.seed!(seed)

    @unpack tlength = time_params

    state_arr = zeros(Float64, size(states, 1), tlength)
    change_arr = zeros(Float64, size(states, 1), tlength)
    jump_arr = zeros(Float64, 9, tlength)

    if retbetaarr == true
        beta_arr = zeros(Float64, tlength)
        seir_wide_mod!(
            state_arr, change_arr, jump_arr, beta_arr, states, dynamics_params,
            time_params;
            type = type, seed = seed,
        )
        return state_arr, change_arr, jump_arr, beta_arr
    else
        seir_wide_mod!(
            state_arr, change_arr, jump_arr, states, dynamics_params,
            time_params;
            type = type,
            seed = seed,
        )
        return state_arr, change_arr, jump_arr
    end
end

"""
    seir_mod!(state_arr, change_arr, jump_arr, states, dynamics_params, trange; tstep, type = "stoch")

The in-place function to run the SEIR model, without producing the transmission rate array.
"""
function seir_wide_mod!(
    state_arr, change_arr, jump_arr, states, dynamics_params, time_params;
    type = "stoch",
    seed = 1234,
)
    Random.seed!(seed)

    for (j, t) in pairs(time_params.trange)
        if j == 1
            state_arr[:, j] = states
            continue
        end

        seir_mod_wide_loop!(
            state_arr, change_arr, jump_arr, j, dynamics_params, time_params, t;
            type = type,
        )
    end

    return nothing
end

"""
    seir_mod!(state_arr, change_arr, jump_arr, beta_arr, states, dynamics_params, trange; tstep, type = "stoch")

The in-place function to run the SEIR model and produce the transmission rate array.
"""
function seir_wide_mod!(
    state_arr, change_arr, jump_arr, beta_arr, states, dynamics_params,
    time_params;
    type = "stoch",
    seed = 1234,
)
    Random.seed!(seed)

    @unpack beta_mean, beta_force = dynamics_params

    for (j, t) in pairs(time_params.trange)
        beta_t = calculate_beta_amp(beta_mean, beta_force, t)
        beta_arr[j] = beta_t

        if j == 1
            state_arr[:, j] = states

            continue
        end

        seir_mod_wide_loop!(
            state_arr, change_arr, jump_arr, j, dynamics_params, time_params, t;
            type = type,
        )
    end

    return nothing
end

"""
    seir_mod_loop!(state_arr, change_arr, jump_arr, j, params, t, tstep; type = type)

The inner loop that is called by `seir_mod!()` function.
"""
function seir_mod_wide_loop!(
    state_arr, change_arr, jump_arr, j, dynamics_params, time_params, t;
    type = type
)
    # Unpack the state variables for easier use
    S = state_arr[1, j - 1]
    E = state_arr[2, j - 1]
    I = state_arr[3, j - 1]
    R = state_arr[4, j - 1]
    N = state_arr[5, j - 1]

    # Unpack the parameters for easier use
    @unpack tstep = time_params
    @unpack beta_mean, beta_force, sigma, gamma, mu, epsilon, R_0 =
        dynamics_params
    beta_t = calculate_beta_amp(beta_mean, beta_force, t)

    # Calculate the rates of each event
    infec_rate = beta_t * S * I        # 1
    latent_rate = sigma * E             # 2
    recov_rate = gamma * I              # 3
    birth_rate = mu * N              # 4
    S_death_rate = mu * S            # 5
    E_death_rate = mu * E            # 6
    I_death_rate = mu * I            # 7
    R_death_rate = mu * R            # 8
    import_rate = epsilon * N / R_0        # 9

    rates = [
        infec_rate, latent_rate, recov_rate, birth_rate, S_death_rate,
        E_death_rate, I_death_rate, R_death_rate, import_rate]

    # Calculate the number of jumps for each event
    if type == "stoch"

        jump_arr[:, j] = map(
            r -> rand(Poisson(r * tstep)),
            rates,
        )
    elseif type == "det"
        jump_arr[:, j] = map(
            r -> r * tstep,
            rates,
        )
    else
        return ("Type must be stoch or det")
    end

    # Calculate the change in each state
    change_arr[1, j] =
        jump_arr[4, j] - jump_arr[1, j] - jump_arr[5, j] - jump_arr[9, j]
    change_arr[2, j] =
        jump_arr[1, j] - jump_arr[2, j] - jump_arr[6, j] + jump_arr[9, j]
    change_arr[3, j] = jump_arr[2, j] - jump_arr[3, j] - jump_arr[7, j]
    change_arr[4, j] = jump_arr[3, j] - jump_arr[8, j]
    change_arr[5, j] = sum(change_arr[1:4, j])

    # Check that the change in each state does not result in a negative state, 
    # and if it is, set the change to the negative of the current state (i.e.
    # set the new state value to 0)
    for state in axes(state_arr, 1)
        if state_arr[state, j - 1] + change_arr[state, j] < 0
            change_arr[state, j] = -state_arr[state, j - 1]
        end
    end

    @. state_arr[:, j] = state_arr[:, j - 1] + change_arr[:, j]

    return nothing
end
































"""
    seir_mod(states, dynamics_params, trange; retbetaarr = false, type = "stoch")

A Tau-leaping SEIR model with commuter style importers.
the model can return an array of transmission rates if `retbetaarr = true`, and can be set to be deterministic or stochastic.

```julia-repl
seir_array, change_array, jump_array, beta_arr = seir_mod(
    init_states, dynamics_params, trange; retbetaarr = true, type = "stoch"
);
```
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
        Float64, time_params.tlength, (size(states, 1) - 1) * 2 + 1
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

    @turbo @. beta_arr = calculate_beta_amp(
        dynamics_params.beta_mean,
        dynamics_params.beta_force,
        time_params.trange,
    )

    for (i, t) in pairs(time_params.trange)
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
    # Calculate the rates of each event
    @views rates[1] = beta_arr[i] .* state_arr[i - 1, 1] .* state_arr[i - 1, 2] # S -> E
    @views rates[2] = dynamics_params.sigma .* state_arr[i - 1, 2] # E -> I
    @views rates[3] = dynamics_params.gamma .* state_arr[i - 1, 3] # I -> R
    @views rates[4] = dynamics_params.mu .* state_arr[i - 1, 5] # Birth -> S
    @views rates[5] = dynamics_params.mu .* state_arr[i - 1, 1] # S -> death
    @views rates[6] = dynamics_params.mu .* state_arr[i - 1, 3] # E -> death
    @views rates[7] = dynamics_params.mu .* state_arr[i - 1, 2] # I -> death
    @views rates[8] = dynamics_params.mu .* state_arr[i - 1, 4] # R -> death
    @views rates[9] =
        dynamics_params.epsilon .* state_arr[i - 1, 5] ./ dynamics_params.R_0 # Import -> S

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
    @views change_arr[i, 4] = jump_arr[i, 3] - jump_arr[i, 8]
    @views change_arr[i, 5] = sum(change_arr[i, 1:4])

    # Check that the change in each state does not result in a negative state,
    # and if it is, set the change to the negative of the current state (i.e.
    # set the new state value to 0)
    for state in axes(state_arr, 2)
        @views if state_arr[i - 1, state] + change_arr[i, state] < 0
            @views change_arr[i, state] = -state_arr[i - 1, state]
        end
    end

    previous_states = @view state_arr[i - 1, :]
    current_changes = @view change_arr[i, :]
    @turbo @. state_arr[i, :] = previous_states + current_changes

    return nothing
end

# end
