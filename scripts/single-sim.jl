#%%
using DrWatson
@quickactivate "OutbreakDetection"

using Random

# include("../src/OutbreakDetection.jl")
using OutbreakDetection


#%%
singlesim_states_p = StateParameters(
    N = 500_000,
    s_prop = 0.1,
    e_prop = 0.01,
    i_prop = 0.01
)

singlesim_time_p = SimTimeParameters(;
    tmin = 0.0, tmax = 365.0 * 100, tstep = 1.0
)

latent_per_days = 8
dur_inf_days = 5
R_0 = 10.0
sigma = 1 / latent_per_days
gamma = 1 / dur_inf_days
mu = 1 / (62.5 * 365)
# beta_mean is the average transmission rate
beta_mean = calculate_beta(R_0, gamma, mu, 1, singlesim_states_p.init_states.N)
# Adjust the scale of the seasonal variation in infectivity i.e. beta_force scales the amplitude of cosine function
beta_force = 0.2
epsilon = calculate_import_rate(mu, R_0, singlesim_states_p.init_states.N)

singlesim_dynamics_p = DynamicsParameters(;
    beta_mean = beta_mean,
    beta_force = beta_force,
    sigma = sigma,
    gamma = gamma,
    mu = mu,
    epsilon = epsilon,
    R_0 = R_0,
)

#%%
seir_array, change_array, jump_array, beta_arr = seir_mod(
    singlesim_states_p.init_states, singlesim_dynamics_p, singlesim_time_p; retbetaarr = true,
    type = "stoch", seed = 1234
);

seir_df = create_sir_df(
    seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N]
)
