#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using JLD2

using OutbreakDetection

#%%
@unpack singlesim_states_p, singlesim_time_p, singlesim_dynamics_p = load("data/singlesim/single-sim_setup.jld2")

#%%
seir_static_mod(
    singlesim_states_p.init_states, singlesim_dynamics_p, singlesim_time_p;
    type = "stoch", seed = 1234,
);
