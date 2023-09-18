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

seir_wide_array, change_wide_array, jump_wide_array, beta_wide_arr = seir_wide_mod(
    singlesim_states_p.init_states, singlesim_dynamics_p, singlesim_time_p;
    retbetaarr = true, type = "stoch", seed = 1234,
);

seir_wide_array' == seir_array

beta_wide_arr == beta_arr

compare_betas = DataFrame("time" => trange, "long" => beta_arr, "wide" => beta_wide_arr)

@subset(compare_betas, :long .!= :wide)

#%%
seir_df = create_sir_df(
    seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N]
)

#%%
jldsave("data/singlesim/single-sim_arrays.jld2"; seir_array, change_array, jump_array, beta_arr, seir_df)
