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

    state_arr = zeros(Int64, time_params.tlength, size(states, 1) + 1)
    beta_arr = zeros(Float64, time_params.tlength)

    seir_mod!(
        state_arr,
        Vector{Int64}(undef, 5),
        Vector{Int64}(undef, 10),
        beta_arr,
        states,
        Vector{Float64}(undef, 6),
        dynamics_params,
        time_params;
        type = type,
        seed = seed,
    )
    return state_arr, beta_arr
end

"""
    seir_mod!(state_arr, change_arr, jump_arr, beta_arr, states, dynamics_params, trange; tstep, type = "stoch")

The in-place function to run the SEIR model and produce the transmission rate array.
"""
function seir_mod!(
    state_arr,
    change_vec,
    jump_vec,
    beta_arr,
    states,
    poisson_rates,
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

    @inbounds for i in eachindex(time_params.trange)
        if i == 1
            state_arr[i, :] .= [states..., 0]
            continue
        end

        seir_mod_loop!(
            state_arr,
            change_vec,
            jump_vec,
            beta_arr,
            poisson_rates,
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
    change_vec,
    jump_vec,
    beta_arr,
    poisson_rates,
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

    @inbounds poisson_rates[1] = beta_t * S * I # Contact: S -> E
    @inbounds poisson_rates[2] = dynamics_params.mu * (1 - dynamics_params.vaccination_coverage) * N # Birth -> S
    @inbounds poisson_rates[3] = dynamics_params.mu * S # S -> death
    @inbounds poisson_rates[4] = dynamics_params.mu * R # R -> death
    @inbounds poisson_rates[5] = dynamics_params.epsilon * N / dynamics_params.R_0 # Import: S -> E
    @inbounds poisson_rates[6] = dynamics_params.mu * dynamics_params.vaccination_coverage * N # Birth -> R

    @simd for r in eachindex(poisson_rates)
        jump_vec[r] = rand(Poisson(poisson_rates[r] * time_params.tstep))
    end

    @inbounds jump_vec[7] = rand(Binomial(E, dynamics_params.sigma * time_params.tstep)) # E -> I
    @inbounds jump_vec[8] = rand(Binomial(E - jump_vec[7], dynamics_params.mu * time_params.tstep)) # E -> death
    @inbounds jump_vec[9] = rand(Binomial(I, dynamics_params.gamma * time_params.tstep)) # I -> R
    @inbounds jump_vec[10] = rand(Binomial(I - jump_vec[9], dynamics_params.mu * time_params.tstep)) # I -> death

    @inbounds contact_inf = jump_vec[1]
    @inbounds S_births = jump_vec[2]
    @inbounds S_death = jump_vec[3]
    @inbounds R_death = jump_vec[4]
    @inbounds import_inf = jump_vec[5]
    @inbounds R_births = jump_vec[6]
    @inbounds latent = jump_vec[7]
    @inbounds E_death = jump_vec[8]
    @inbounds recovery = jump_vec[9]
    @inbounds I_death = jump_vec[10]

    # Calculate the change in each state
    @inbounds change_vec[1] = S_births - (contact_inf + import_inf + S_death)
    @inbounds change_vec[2] = (contact_inf + import_inf) - (latent + E_death)
    @inbounds change_vec[3] = latent - (recovery + I_death)
    @inbounds change_vec[4] = (recovery + R_births) - R_death
    @inbounds change_vec[5] = sum(@view(change_vec[1:4]))

    @simd for j in 1:5
        @views state_arr[i, j] = state_arr[i - 1, j] + change_vec[j]
    end
    @inbounds state_arr[i, end] = jump_vec[1]

    return nothing
end

# end
