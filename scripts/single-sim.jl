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

#%%
seir_vec, beta_vec = seir_mod(
    init_states_static, singlesim_dynamics_p, singlesim_time_p;
    seed = 1234,
)

seir_array = convert_svec_to_arr(seir_vec; reinterpret_dims = (6, length(singlesim_time_p.trange)), reorder_inds = (2, 1))

seir_df = create_sir_df(
    seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N, :incidence]
)

#%%
jldsave(
    "data/singlesim/single-sim_arrays.jld2";
    seir_array,
    change_array,
    jump_array,
    beta_arr,
    seir_df,
)
