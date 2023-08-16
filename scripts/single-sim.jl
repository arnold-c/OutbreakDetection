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
N = 5e5
s = 0.1
e = 0.01
i = 0.01
r = 1.0 - (s + e + i)
u₀ = convert.(Int64, [N * s, N * e, N * i, N * r, N])
singlesim_time_p = SimTimeParameters(0.0, 365.0 * 100, 1.0)

latent_per = 8
dur_inf = 5
R₀ = 10.0
sigma = 1 / latent_per
gamma = 1 / dur_inf
mu = 1 / (62.5 * 365)
# beta_mean is the average transmission rate
beta_mean = calculate_beta(R₀, gamma, mu, 1, N)
# Adjust the scale of the seasonal variation in infectivity i.e. beta_force scales the amplitude of cosine function
beta_force = 0.2
ε = (1.06 * mu * (R₀ - 1)) / sqrt(N) # Commuter imports - see p210 Keeling & Rohani
p = (beta_mean, beta_force, sigma, gamma, mu, ε, R₀)

Random.seed!(1234)

#%%
seir_array, change_array, jump_array, beta_arr = seir_mod(
    u₀, p, trange; retbetaarr = true, type = "stoch"
);

seir_df = create_sir_df(
    seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N]
)
