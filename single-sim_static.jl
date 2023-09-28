#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using JLD2
using StaticArrays

using OutbreakDetection

#%%
@unpack singlesim_states_p, singlesim_time_p, singlesim_dynamics_p = load("data/singlesim/single-sim_setup.jld2")

static_states = @SVector [Float64(singlesim_states_p.init_states[i]) for i in 1:5]
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

reinterpret(Float64, static_seir)

@chain reshape(reinterpret(Float64, static_seir), (5, length(singlesim_time_p.trange)))' begin
    Array(_)
    create_sir_df(_, singlesim_time_p.trange, [:S, :E, :I, :R, :N])
    draw_sir_plot(_, labels = ["S", "E", "I", "R", "N"])
end

#%%
seir_arr = seir_mod(singlesim_states_p.init_states, singlesim_dynamics_p, singlesim_time_p; type = "stoch", seed = 1234)[1]

static_seir == seir_arr

seir_arr[end, :]
static_seir[end]

#%%
@benchmark seir_mod(singlesim_states_p.init_states, singlesim_dynamics_p, singlesim_time_p; type = "stoch", seed = 1234)

@benchmark seir_static_mod(
    static_states, singlesim_dynamics_p, singlesim_time_p;
    type = "stoch", seed = 1234,
)
