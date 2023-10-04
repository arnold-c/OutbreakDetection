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
seir_array, beta_arr = seir_mod(
    singlesim_states_p.init_states, singlesim_dynamics_p, singlesim_time_p;
    type = "stoch", seed = 1234,
);

seir_df = create_sir_df(
    seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N, :incidence]
)

#%%
@benchmark seir_mod(
    $singlesim_states_p.init_states, $singlesim_dynamics_p, $singlesim_time_p;
    type = "stoch", seed = $1234,
)

@benchmark seir_mod!(
    $seir_array,
    # $Vector{Int64}(undef, 5),
    # $Vector{Int64}(undef, 10),
    $MVector{5,Int64}(undef),
    $MVector{10,Int64}(undef),
    $beta_arr,
    $singlesim_states_p.init_states,
    # $Vector{Float64}(undef, 6),
    $MVector{6,Float64}(undef),
    $singlesim_dynamics_p,
    $singlesim_time_p;
    type = "stoch", seed = $1234,
)

init_states_static = @SVector [
    singlesim_states_p.init_states[1],
    singlesim_states_p.init_states[2],
    singlesim_states_p.init_states[3],
    singlesim_states_p.init_states[4],
    singlesim_states_p.init_states[5],
    0
]
seir_vec = Vector{typeof(init_states_static)}(
    undef, length(singlesim_time_p.trange)
)

seir_mod_static!(
    seir_vec,
    MVector{5,Int64}(undef),
    MVector{10,Int64}(undef),
    beta_arr,
    init_states_static,
    MVector{6,Float64}(undef),
    singlesim_dynamics_p,
    singlesim_time_p;
    type = "stoch", seed = 1234,
)

@benchmark seir_mod_static!(
    $seir_vec,
    $MVector{5,Int64}(undef),
    $MVector{10,Int64}(undef),
    $beta_arr,
    $init_states_static,
    $MVector{6,Float64}(undef),
    $singlesim_dynamics_p,
    $singlesim_time_p;
    type = "stoch", seed = $1234,
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
