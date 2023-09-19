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
    # Create views of the state variables for easier use
    @views S = state_arr[i - 1, 1]
    @views E = state_arr[i - 1, 2]
    @views I = state_arr[i - 1, 3]
    @views R = state_arr[i - 1, 4]
    @views N = state_arr[i - 1, 5]
    @views beta_t = beta_arr[i]

    # Calculate the rates of each event
    rates[1] = beta_t * S * I              # S -> E
    rates[2] = dynamics_params.sigma * E        # E -> I
    rates[3] = dynamics_params.gamma * I        # I -> R
    rates[4] = dynamics_params.mu * N           # Birth -> S
    rates[5] = dynamics_params.mu * S           # S -> death
    rates[6] = dynamics_params.mu * E           # E -> death
    rates[7] = dynamics_params.mu * I           # I -> death
    rates[8] = dynamics_params.mu * R           # R -> death
    rates[9] = dynamics_params.epsilon * N / dynamics_params.R_0        # Import -> S

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
    @views change_arr[i, 1] = jump_arr[i, 4] - jump_arr[i, 1] - jump_arr[i, 5] - jump_arr[i, 9]
    @views change_arr[i, 2] = jump_arr[i, 1] - jump_arr[i, 2] - jump_arr[i, 6] + jump_arr[i, 9]
    @views change_arr[i, 3] = jump_arr[i, 2] - jump_arr[i, 3] - jump_arr[i, 7]
    @views change_arr[i, 4] = jump_arr[i, 3] - jump_arr[i, 8]

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
    @turbo @. state_arr[i, :] = previous_states + current_changes

    return nothing
end

# end
