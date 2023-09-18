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
    type = "det", seed = 1234,
);
seir_wide_array' == seir_array

seir_wide_array, change_wide_array, jump_wide_array, beta_wide_arr = seir_wide_mod(
    singlesim_states_p.init_states, singlesim_dynamics_p, singlesim_time_p;
    retbetaarr = true, type = "det", seed = 1234,
);


change_array[1:10, :]
change_wide_array'[1:10, :]


compare_jumps = DataFrame("time" => trange, "long" => jump_array[:, 1], "wide" => jump_wide_array[1, :])

@subset(compare_jumps, :long .!= :wide)

#%%
seir_df = create_sir_df(
    seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N]
)

#%%
jldsave("data/singlesim/single-sim_arrays.jld2"; seir_array, change_array, jump_array, beta_arr, seir_df)
