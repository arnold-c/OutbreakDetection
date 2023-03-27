"""
This is a stochastic differential equation SIR model using the ModelingToolkit interface.
The noise is diagonal (i.e., every function in the system gets a different random number) and the noise is added to the S, I and R variables.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using DifferentialEquations, ModelingToolkit, Statistics
using DataFrames, DataFramesMeta, LinearAlgebra, Symbolics
using CairoMakie, AlgebraOfGraphics

CairoMakie.activate!()
set_aog_theme!()

#%%
includet(srcdir("Julia/transmission-functions.jl")) # Revise will keep track of file changes and reload functions as necessary
includet(srcdir("Julia/plotting-functions.jl"))
includet(srcdir("Julia/cleaning-functions.jl"))

#%%
δt = 0.005
tlower = 0.0
tmax = 250.0
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

#%%
@parameters β γ μ
@variables t S(t) I(t) R(t)
D = Differential(t)

#%%
nic = 1 # number of infective compartments
nac = 1 # number of age classes

#%%
eqs = [D(S) ~ -β * S * I + μ * (I + R), D(I) ~ β * S * I - γ * I - μ * I,
    D(R) ~ γ * I - μ * R]

# Noise is sampled from a normal distribution with mean 0 and standard deviation of δt time-step, and multiplied by the value specified in the noise equation, e.g., 0.5 + 0.01 * S * 𝒩(0, δt)
noise_eqs = [
    0.5 + 0.01 * S,
    0.5 + 0.01 * I,
    0.5 + 0.01 * R,
]

#%%
@named ode = ODESystem(eqs, t, [S, I, R], [β, γ, μ]; tspan = tspan)
ode_simple = structural_simplify(ode)

@named sde = SDESystem(ode, noise_eqs)

#%%
R₀ = 10.0
p = Dict(γ => 1 / 8, μ => 1 / (62 * 365))
u₀ = Dict(S => 999.0, I => 1.0, R => 0.0)
push!(
    p,
    β => calculate_beta(ode_simple, 1, 1, R₀, p, [1.0], [sum(values(u₀))]),
)

#%%
sde_prob = SDEProblem(sde, u₀, tspan, p)
sde_sol = solve(sde_prob, SOSRI(); saveat = δt)
sde_sol_df = create_sir_df(sde_sol)

#%%
colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(sde_sol_df; colors = colors)

#%%
nsims = 1000

sir_array = zeros(5, tlength)
all_sims_array = fill(NaN, 4, tlength, nsims)

quantiles = [0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975]
sim_quantiles = zeros(Float64, length(quantiles), tlength, 4)

#%%
ensemble_prob = EnsembleProblem(sde_prob)
ensemble_sol = solve(ensemble_prob, EnsembleThreads(); trajectories = nsims, saveat = δt)

create_sir_all_sims_array!(ensemble_sol, nsims)

create_sir_all_sim_quantiles!(; quantiles = quantiles)

create_sir_quantiles_plot(; lower = 0.025, upper = 0.975, quantiles = quantiles)