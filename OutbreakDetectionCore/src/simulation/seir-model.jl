# module SEIRModel
#
# export calculate_beta_amp, seir_mod, seir_mod!, seir_mod_loop!
#
"""
This is a simulation of an SIR model that uses Tau-leaping, with commuter
imports. All jumps are manually defined.
"""

export seir_mod,
    seir_mod!

"""
    seir_mod(states, dynamics_params, trange; tstep, type = "stoch")

The in-place function to run the SEIR model with a vaccinations going directly to the R compartment and produce the transmission rate array.
"""
function seir_mod(
        states::StaticArrays.SVector{5, Int64},
        dynamics_params::DynamicsParameters,
        time_params::SimTimeParameters;
        seed::Int64 = 1234
    )
    Random.seed!(seed)

    tlength = time_params.tlength
    state_vec = Vector{typeof(states)}(undef, tlength)
    inc_vec = Vector{Int64}(undef, tlength)
    Reff_vec = Vector{Float64}(undef, tlength)
    beta_vec = Vector{Float64}(undef, tlength)

    calculate_beta_amp!(
        beta_vec,
        dynamics_params,
        time_params
    )

    seir_mod!(
        state_vec,
        inc_vec,
        Reff_vec,
        beta_vec,
        states,
        dynamics_params,
        time_params,
    )

    return SEIRRun(
        state_vec,
        inc_vec,
        Reff_vec,
    )
end

function seir_mod(
        states::StaticArrays.SVector{5, Int64},
        dynamics_params::DynamicsParameters,
        beta_vec::Vector{Float64},
        time_params::SimTimeParameters;
        seed::Int64 = 1234
    )
    Random.seed!(seed)

    tlength = time_params.tlength
    state_vec = Vector{typeof(states)}(undef, tlength)
    inc_vec = Vector{Int64}(undef, tlength)
    Reff_vec = Vector{Float64}(undef, tlength)

    seir_mod!(
        state_vec,
        inc_vec,
        Reff_vec,
        beta_vec,
        states,
        dynamics_params,
        time_params,
    )

    return SEIRRun(
        state_vec,
        inc_vec,
        Reff_vec,
    )
end


"""
    seir_mod!(state_arr, change_arr, jump_arr, beta_arr, states, dynamics_params, trange; tstep, type = "stoch")

The in-place function to run the SEIR model and produce the transmission rate array.
"""
function seir_mod!(
        state_vec::ASV,
        inc_vec::AI,
        Reff_vec::AF1,
        beta_vec::AF2,
        states::StaticArrays.SVector{5, Int64},
        dynamics_params::DynamicsParameters,
        time_params::SimTimeParameters,
    ) where {
        ASV <: AbstractVector{StaticArrays.SVector{5, Int64}},
        AI <: AbstractVector{Int64},
        AF1 <: AbstractVector{Float64},
        AF2 <: AbstractVector{Float64},
    }

    @inbounds begin
        mu = dynamics_params.mu
        epsilon = dynamics_params.epsilon
        sigma = dynamics_params.sigma
        gamma = dynamics_params.gamma
        timestep = time_params.tstep
        R_0 = dynamics_params.R_0
        tlength = length(inc_vec)

        state_vec[1] = states
        inc_vec[1] = 0
        Reff_vec[1] = calculateReffective(
            beta_vec[1],
            dynamics_params,
            states[1],
            states[5],
        )
    end

    # Pre-compute constant values to avoid repeated calculations
    mu_timestep = mu * timestep
    sigma_timestep = sigma * timestep
    gamma_timestep = gamma * timestep
    epsilon_timestep = epsilon * timestep

    @inbounds for i in 2:(tlength)
        vaccination_coverage = dynamics_params.vaccination_coverage

        state_vec[i], inc_vec[i] = seir_mod_loop(
            state_vec[i - 1],
            beta_vec[i],
            mu_timestep,
            epsilon_timestep,
            sigma_timestep,
            gamma_timestep,
            R_0,
            vaccination_coverage,
            timestep,
        )
        Reff_vec[i] = calculateReffective(
            beta_vec[i],
            dynamics_params,
            state_vec[i][1],
            state_vec[i][5],
        )
    end

    return nothing
end

"""
    seir_mod_loop!(state_arr, change_arr, jump_arr, j, params, t, tstep; type = type)

The inner loop that is called by `seir_mod!()` function.
"""
function seir_mod_loop(
        state_vec::StaticArrays.SVector{5, Int64},
        beta_t::Float64,
        mu_timestep::Float64,
        epsilon_timestep::Float64,
        sigma_timestep::Float64,
        gamma_timestep::Float64,
        R_0::Float64,
        vaccination_coverage::Float64,
        timestep::Float64,
    )::Tuple{StaticArrays.SVector{5, Int64}, Int64}

    @inbounds begin
        S = state_vec[1]
        E = state_vec[2]
        I = state_vec[3]
        R = state_vec[4]
        N = state_vec[5]

        contact_inf = Random.rand(
            Distributions.Poisson(beta_t * S * I * timestep)
        ) # Contact: S -> E
        S_births = Random.rand(
            Distributions.Poisson(
                (1 - vaccination_coverage) * N * mu_timestep
            ),
        ) # Birth -> S
        S_death = Random.rand(Distributions.Poisson(S * mu_timestep)) # S -> death
        R_death = Random.rand(Distributions.Poisson(R * mu_timestep)) # R -> death
        import_inf = Random.rand(
            Distributions.Poisson(epsilon_timestep * N / R_0)
        ) # Import: S -> E
        R_births = Random.rand(
            Distributions.Poisson(vaccination_coverage * N * mu_timestep)
        ) # Birth -> R
        latent = Random.rand(Distributions.Binomial(E, sigma_timestep)) # E -> I
        E_death = Random.rand(Distributions.Binomial(E - latent, mu_timestep)) # E -> death
        recovery = Random.rand(Distributions.Binomial(I, gamma_timestep)) # I -> R
        I_death = Random.rand(Distributions.Binomial(I - recovery, mu_timestep)) # I -> death

        # Calculate changes
        dS = S_births - (contact_inf + import_inf + S_death)
        dE = (contact_inf + import_inf) - (latent + E_death)
        dI = latent - (recovery + I_death)
        dR = (recovery + R_births) - R_death
        dN = dS + dE + dI + dR
    end

    return (
        StaticArrays.SVector(S + dS, E + dE, I + dI, R + dR, N + dN),
        contact_inf,
    )
end

"""
    convert_rate_to_prob(rate)

A simple function to ensure that no rate is greater than 1.
Rests on the assumption that the rate and timestep are both sufficiently small,
therefore, the probability ≈ rate within a very small tolerance, so just return
the rate.
"""
@inline function convert_rate_to_prob(rate)
    return min(rate, 1.0)
end

@inline function _smart_transition(n::Int64, prob::Float64)
    rate = n * prob
    return _smart_transition(n, prob, rate)
end

# Smart transition function that chooses between Poisson and Binomial
@inline function _smart_transition(n::Int64, prob::Float64, rate::Float64)
    # Use Poisson when:
    # - Population is large (n > 30)
    # - Expected value is small (rate < 10)
    # This gives good approximation with better performance

    if n == 0
        return 0
    end

    if n > 30 && rate < 10.0
        # Use Poisson approximation for large n, small rate
        return min(rand(Distributions.Poisson(rate)), n)  # Cap at population size
    else
        # Use exact Binomial for small n or large rate
        return rand(Distributions.Binomial(n, prob))
    end
end

function convert_svec_to_matrix(svec)
    arr = Matrix{Int64}(undef, size(svec, 1), length(svec[1]))
    convert_svec_to_matrix!(arr, svec)
    return arr
end
function convert_svec_to_matrix!(arr, svec)
    @inbounds for state in eachindex(svec[1]), time in axes(svec, 1)
        arr[time, state] = svec[time][state]
    end
    return nothing
end

function convert_svec_to_array(svec)
    arr = Array{Int64}(undef, size(svec, 1), length(svec[1]), size(svec, 2))
    @inbounds for sim in axes(svec, 2)
        convert_svec_to_matrix!(@view(arr[:, :, sim]), @view(svec[:, sim]))
    end
    return arr
end
