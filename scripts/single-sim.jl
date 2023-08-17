#%%
using DrWatson
@quickactivate "OutbreakDetection"

using Random

#%%
# Revise will keep track of file change_arr and reload functions as necessary
includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(funsdir("transmission-functions.jl"))
includet(funsdir("plotting-functions.jl"))
includet(funsdir("cleaning-functions.jl"))
includet(funsdir("SEIR-model.jl"))
includet(funsdir("structs.jl"))

#%%
N = 500_000
s_init_prop = 0.1
e_init_prop = 0.01
i_init_prop = 0.01
r_init_prop = 1.0 - (s_init_prop + e_init_prop + i)
init_states =
    convert.(
        Int64,
        [N * s_init_prop, N * e_init_prop, N * i_init_prop, N * r_init_prop, N],
    )
singlesim_time_p = SimTimeParameters(0.0, 365.0 * 100, 1.0)

latent_per = 8
dur_inf = 5
R_0 = 10.0
sigma = 1 / latent_per
gamma = 1 / dur_inf
mu = 1 / (62.5 * 365)
# beta_mean is the average transmission rate
beta_mean = calculate_beta(R_0, gamma, mu, 1, N)
# Adjust the scale of the seasonal variation in infectivity i.e. beta_force scales the amplitude of cosine function
beta_force = 0.2
epsilon = (1.06 * mu * (R_0 - 1)) / sqrt(N) # Commuter imports - see p210 Keeling & Rohani
singlesim_p = (beta_mean, beta_force, sigma, gamma, mu, epsilon, R_0)

Random.seed!(1234)

#%%
seir_array, change_array, jump_array, beta_arr = seir_mod(
    init_states, singlesim_p, trange; retbetaarr = true, type = "stoch"
);

seir_df = create_sir_df(
    seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N]
)
