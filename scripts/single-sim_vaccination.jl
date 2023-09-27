#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using JLD2

using OutbreakDetection

include(srcdir("makie-plotting-setup.jl"))

#%%
@unpack singlesim_states_p, singlesim_time_p, singlesim_dynamics_p = load("data/singlesim/single-sim_setup.jld2")
@unpack seir_array = load("data/singlesim/single-sim_arrays.jld2")

#%%
seirv_array, changev_array, jumpv_array, betav_arr = seirv_mod(
    singlesim_states_p.init_states,
    DynamicsParameters(
        singlesim_dynamics_p.beta_mean,
        singlesim_dynamics_p.beta_force,
        singlesim_dynamics_p.sigma,
        singlesim_dynamics_p.gamma,
        singlesim_dynamics_p.mu,
        singlesim_dynamics_p.annual_births_per_k,
        singlesim_dynamics_p.epsilon,
        singlesim_dynamics_p.R_0,
        0.0
    ),
    singlesim_time_p; type = "stoch",
    seed = 1234,
);

seir_array == seirv_array

#%%
seirv_df = create_sir_df(
    seirv_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N]
)

#%%
singlesim_vaccination_timeseries_plot = draw_sir_plot(
    seirv_df;
    annual = true,
    colors = seircolors,
    labels = seir_state_labels
)

#%%
draw_sir_plot(
    create_sir_df(seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N]);
    annual = true,
    colors = seircolors,
    labels = seir_state_labels
)
