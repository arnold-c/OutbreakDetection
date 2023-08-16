"""
This is a simulation of an SIR model that uses Tau-leaping, with commuter
imports. All jumps are manually defined.
"""

using DrWatson
@quickactivate "OutbreakDetection"

using DifferentialEquations
using Statistics
using Distributions
using Random

"""
    calculate_beta_amp(beta_mean, beta_force, t)

Calculate the amplitude of the transmission rate beta as a function of time.
`beta_mean` is the mean transmission rate, `beta_force` is the amplitude of the cosine function.
"""
function calculate_beta_amp(beta_mean, beta_force, t)
    return beta_mean * (1 + beta_force * cos(2pi * t / 365))
end

"""
    seir_mod(u, p, trange; retbetaarr = false, type = "stoch")

A Tau-leaping SEIR model with commuter style importers.
the model can return an array of transmission rates if `retbetaarr = true`, and can be set to be deterministic or stochastic.

```julia-repl
seir_array, change_array, jump_array, beta_arr = seir_mod(
    init_states, p, trange; retbetaarr = true, type = "stoch"
);
```
"""
function seir_mod(u, p, trange; retbetaarr = false, type = "stoch", seed = 1234)
    tlength = length(trange)
    tstep = step(trange)

    state_arr = zeros(Float64, size(u, 1), tlength)

    change_arr = zeros(Float64, size(u, 1), tlength)

    jump_arr = zeros(Float64, 9, tlength)

    if retbetaarr == true
        beta_arr = zeros(Float64, tlength)
        seir_mod!(
            state_arr, change_arr, jump_arr, beta_arr, u, p, trange; tstep = tstep,
            type = type, seed = seed,
        )
        return state_arr, change_arr, jump_arr, beta_arr
    else
        seir_mod!(
            state_arr, change_arr, jump_arr, u, p, trange; tstep = tstep, type = type,
            seed = seed,
        )
        return state_arr, change_arr, jump_arr
    end
end

"""
    seir_mod!(state_arr, change_arr, jump_arr, u, p, trange; tstep, type = "stoch")

The in-place function to run the SEIR model, without producing the transmission rate array.
"""
function seir_mod!(
    state_arr, change_arr, jump_arr, u, p, trange; tstep, type = "stoch",
    seed = 1234,
)
    for (j, t) in pairs(trange)
        if j == 1
            state_arr[:, j] = u
            continue
        end

        seir_mod_loop!(
            state_arr, change_arr, jump_arr, j, p, t, dt; type = type,
            seed = seed,
        )
    end

    return nothing
end

"""
    seir_mod!(state_arr, change_arr, jump_arr, beta_arr, u, p, trange; tstep, type = "stoch")

The in-place function to run the SEIR model and produce the transmission rate array.
"""
function seir_mod!(
    state_arr, change_arr, jump_arr, beta_arr, u, p, trange; tstep, type = "stoch",
    seed = 1234,
)
    beta_mean, beta_force = p

    for (j, t) in pairs(trange)
        beta_t = calculate_beta_amp(beta_mean, beta_force, t)
        beta_arr[j] = beta_t

        if j == 1
            state_arr[:, j] = u

            continue
        end

        seir_mod_loop!(
            state_arr, change_arr, jump_arr, j, p, t, dt; type = type,
            seed = seed,
        )
    end

    return nothing
end

"""
    seir_mod_loop!(state_arr, change_arr, jump_arr, j, p, t, dt; type = type)

The inner loop that is called by `seir_mod!()` function.
"""
function seir_mod_loop!(
    state_arr, change_arr, jump_arr, j, p, t, dt; type = type, seed = 1234
)
    # Unpack the state variables for easier use
    S = state_arr[1, j - 1]
    E = state_arr[2, j - 1]
    I = state_arr[3, j - 1]
    R = state_arr[4, j - 1]
    N = state_arr[5, j - 1]

    # Unpack the parameters for easier use
    beta_mean, beta_force, sigma, gamma, mu, epsilon, R_0 = p
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
        Random.seed!(seed)

        jump_arr[:, j] = map(
            r -> rand(Poisson(r * dt)),
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
    # and if it is, set the change to the negative of the current state
    for state in axes(state_arr, 1)
        if state_arr[state, j - 1] + change_arr[state, j] < 0
            change_arr[state, j] = -state_arr[state, j - 1]
        end
    end

    @. state_arr[:, j] = state_arr[:, j - 1] + change_arr[:, j]

    return nothing
end
