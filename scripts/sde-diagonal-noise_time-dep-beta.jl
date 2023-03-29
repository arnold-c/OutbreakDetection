"""
This is a stochastic differential equation SIR model using the ModelingToolkit interface.
The noise is diagonal (i.e., every function in the system gets a different random number) and the noise is added to the S, I and R variables.
A log transformation is applied to β to ensure that it is positive when random walks are incorporated.
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
δt = 0.01
tlower = 0.0
tmax = 250.0
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

#%%
function sir_mod!(du, u, p, t)
    S, I, R, log_β = u
    γ, μ, σ = p

    β = exp(log_β)

    du[1] = -β * S * I + μ * (I + R)
    du[2] = β * S * I - γ * I - μ * I
    du[3] = γ * I - μ * R
    du[4] = 0

    return nothing
end

function noise!(du, u, p, t)
    S, I, R, log_β = u
    γ, μ, σ = p

    du[1] = 0.0
    du[2] = 0.0
    du[3] = 0.0
    du[4] = σ

    return nothing
end

#%%
u₀ = [999.0, 1.0, 0.0]

R₀ = 10.0
γ = 1 / 8
μ = 1 / (62 * 365)
σ = 0.05
β = calculate_beta(R₀, γ, μ, ones(1, 1), [sum(u₀)])

p = [γ, μ, σ]

push!(u₀, log(β))

#%%
sde_prob = SDEProblem(sir_mod!, noise!, u₀, tspan, p)
sde_sol = solve(sde_prob, SOSRI(); saveat = δt)
sde_sol_df, sde_sol_beta = create_sir_df(sde_sol; β = true)
sde_sol_beta.beta = exp.(sde_sol_beta.beta)

#%%
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Time (days)", ylabel = "β")

lines!(ax, sde_sol_beta.time, sde_sol_beta.beta; linewidth = 2)

fig

#%%
colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple", "green"]

create_sir_plot(sde_sol_df; colors = colors)

#%%
nsims = 1000

sir_array = zeros(5, tlength)
all_sims_array = fill(NaN, 5, tlength, nsims)

quantiles = [0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975]

#%%
ensemble_prob = EnsembleProblem(sde_prob)
ensemble_sol = solve(ensemble_prob, EnsembleThreads(); trajectories = nsims, saveat = δt)

create_sir_all_sims_array!(ensemble_sol, nsims; β = true)

all_sims_array[5, :, :] = exp.(all_sims_array[5, :, :])

sim_quantiles = zeros(Float64, length(quantiles), tlength, 5)
create_sir_all_sim_quantiles!(; quantiles = quantiles)

create_sir_quantiles_plot(; lower = 0.025, upper = 0.975, quantiles = quantiles)

create_beta_quantiles_plot(; lower = 0.025, upper = 0.975, quantiles = quantiles)