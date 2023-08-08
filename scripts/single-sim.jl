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

#%%
N = 5e5
s = 0.1
e = 0.01
i = 0.01
r = 1.0 - (s + e + i)
u₀ = convert.(Int64, [N * s, N * e, N * i, N * r, N])
τ = 1.0
tlower = 0.0
tmax = 365.0 * 100
trange = tlower:τ:tmax
tlength = length(trange)

latent_per = 8
dur_inf = 5
R₀ = 10.0
σ = 1 / latent_per
γ = 1 / dur_inf
μ = 1 / (62.5 * 365)
# β₀ is the average transmission rate
β₀ = calculate_beta(R₀, γ, μ, 1, N)
# Adjust the scale of the seasonal variation in infectivity i.e. β₁ scales the amplitude of cosine function
β₁ = 0.2
ε = (1.06 * μ * (R₀ - 1)) / sqrt(N) # Commuter imports - see p210 Keeling & Rohani
p = (β₀, β₁, σ, γ, μ, ε, R₀)

Random.seed!(1234)

#%%
seir_array, change_array, jump_array, β_arr = seir_mod(
    u₀, p, trange; retβarr = true, type = "stoch"
);

seir_df = create_sir_df(seir_array, trange, [:S, :E, :I, :R, :N])

seircolors = ["dodgerblue4", "green", "firebrick3", "chocolate2", "purple"]
state_labels = ["S", "E", "I", "R", "N"]

