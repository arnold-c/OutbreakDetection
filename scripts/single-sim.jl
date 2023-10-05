#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using JLD2

using OutbreakDetection

#%%
@unpack singlesim_states_p, singlesim_time_p, singlesim_dynamics_p = load(
    "data/singlesim/single-sim_setup.jld2"
)

#%%
seir_vec, inc_vec, beta_vec = seir_mod(
    singlesim_states_p.init_states,
    singlesim_dynamics_p,
    singlesim_time_p; seed = 1234,
)

seir_array = convert_svec_to_arr(seir_vec; reinterpret_dims = (5, length(singlesim_time_p.trange)), reorder_inds = (2, 1))

seir_df = create_sir_df(
    seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N]
)

#%%
jldsave(
    "data/singlesim/single-sim_arrays.jld2";
    seir_array,
    inc_vec,
    beta_vec,
    seir_df,
)
