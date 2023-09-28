#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using JLD2
using StaticArrays

using OutbreakDetection

#%%
@unpack singlesim_states_p, singlesim_time_p, singlesim_dynamics_p = load(
    "data/singlesim/single-sim_setup.jld2"
)

static_states = @SVector [
    Float64(singlesim_states_p.init_states[i]) for i in 1:5
]
typeof(static_states)

#%%
static_seir, static_change, static_jump, static_beta = seir_static_mod(
    static_states, singlesim_dynamics_p, singlesim_time_p;
    type = "stoch", seed = 1234,
);
static_seir

static_change
static_jump

lines(static_beta[1:720])

draw_sir_plot(
    create_sir_df(static_seir, singlesim_time_p.trange, [:S, :E, :I, :R, :N]);
    labels = ["S", "E", "I", "R", "N"],
)

# @chain reshape(reinterpret(Float64, static_seir), (5, length(singlesim_time_p.trange))) begin
#     permutedims(_, (2, 1))
# Array(_)
# create_sir_df(_, singlesim_time_p.trange, [:S, :E, :I, :R, :N])
# draw_sir_plot(_, labels = ["S", "E", "I", "R", "N"])
# end

#%%
seir_arr = seir_mod(
    singlesim_states_p.init_states,
    singlesim_dynamics_p,
    singlesim_time_p;
    type = "stoch",
    seed = 1234,
)[1]

static_seir == seir_arr

seir_arr[end, :]
static_seir[end]

#%%
@benchmark seir_mod(
    singlesim_states_p.init_states,
    singlesim_dynamics_p,
    singlesim_time_p;
    type = "stoch",
    seed = 1234,
)

@benchmark seir_static_mod(
    static_states, singlesim_dynamics_p, singlesim_time_p;
    type = "stoch", seed = 1234,
)

#%%
ProfileCanvas.@profview seir_mod(
    singlesim_states_p.init_states,
    singlesim_dynamics_p,
    singlesim_time_p;
    type = "stoch",
    seed = 1234,
)

ProfileCanvas.@profview seir_static_mod(
    static_states, singlesim_dynamics_p, singlesim_time_p;
    type = "stoch", seed = 1234,
)

#%%
function test_regular_ensemble(states, dynamics_params, time_params, nsims)
    seir_arr = Array{Float64}(
        undef, length(time_params.trange), length(states), nsims
    )
    change_arr = similar(seir_arr)
    jump_arr = Array{Float64}(undef, length(time_params.trange), 10, nsims)
    beta_arr = Vector{Float64}(undef, length(time_params.trange))

    for sim in 1:nsims
        seir = @view seir_arr[:, :, sim]
        change = @view change_arr[:, :, sim]
        jump = @view jump_arr[:, :, sim]

        seir_mod!(
            seir,
            change,
            jump,
            beta_arr,
            states,
            Vector{Float64}(undef, 10),
            dynamics_params,
            time_params;
            type = "stoch",
            seed = 1234 + sim,
        )
    end

    return seir_arr
end

@benchmark test_regular_ensemble(
    $singlesim_states_p.init_states, $singlesim_dynamics_p, $singlesim_time_p,
    $1000,
)

#%%
function test_static_ensemble(states, dynamics_params, time_params, nsims)
    seir_arr = Array{typeof(states)}(undef, length(time_params.trange), nsims)
    change_arr = similar(seir_arr)
    jump_arr = Array{SVector{10,eltype(states)}}(
        undef, time_params.tlength, nsims
    )
    @inbounds beta_vec = map(
        t -> calculate_beta_amp(
            dynamics_params.beta_mean, dynamics_params.beta_force, t
        ),
        time_params.trange,
    )

    @simd for sim in 1:nsims
        Random.seed!(1234 + sim)

        @inbounds seir_arr[1, sim] = states

        @inbounds for i in 1:(length(time_params.trange) - 1)
            seir_arr[i + 1, sim], change_arr[i + 1, sim], jump_arr[i + 1, sim] = seir_static_mod_loop(
                beta_vec[i],
                Vector{eltype(states)}(undef, length(states) - 1),
                seir_arr[i, sim],
                dynamics_params,
                time_params;
                type = "stoch",
            )
        end
    end

    return permutedims(reshape(reinterpret(Float64, seir_arr), (5, length(singlesim_time_p.trange), nsims)), (2, 1, 3))
end

@benchmark test_static_ensemble($static_states, $singlesim_dynamics_p, $singlesim_time_p, $1000)
