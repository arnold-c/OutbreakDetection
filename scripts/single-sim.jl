#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using JLD2

using OutbreakDetection

#%%
@unpack singlesim_states_p, singlesim_time_p, singlesim_dynamics_p = load("data/singlesim/single-sim_setup.jld2")

#%%
seir_array, change_array, jump_array, beta_arr = seir_mod(
    singlesim_states_p.init_states, singlesim_dynamics_p, singlesim_time_p;
    type = "stoch", seed = 1234,
);

seir_df = create_sir_df(
    seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N]
)

#%%
@benchmark seir_mod(
    $singlesim_states_p.init_states, $singlesim_dynamics_p, $singlesim_time_p;
    type = "stoch", seed = $1234,
)

@benchmark seir_mod!(
    $seir_array,
    $change_array,
    $jump_array,
    $beta_arr,
    $singlesim_states_p.init_states,
    $Vector{Float64}(undef, 6),
    $singlesim_dynamics_p,
    $singlesim_time_p;
    type = "stoch", seed = $1234,
)



#%%
jldsave("data/singlesim/single-sim_arrays.jld2"; seir_array, change_array, jump_array, beta_arr, seir_df)
