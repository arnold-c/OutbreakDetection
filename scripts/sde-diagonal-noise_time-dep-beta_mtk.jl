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
δt = 0.01
tlower = 0.0
tmax = 250.0
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

#%%
@parameters γ μ σ
@variables t S(t) I(t) R(t) β(t)
D = Differential(t)

#%%
nic = 1 # number of infective compartments
nac = 1 # number of age classes

#%%
eqs = [
    D(S) ~ -exp(β) * S * I + μ * (I + R),
    D(I) ~ exp(β) * S * I - γ * I - μ * I,
    D(R) ~ γ * I - μ * R,
]

trans_eqs = [D(β) ~ 0.0]

noise_eqs = [
    0.0,
    0.0,
    0.0,
    σ,
]

#%%
@named ode = ODESystem([eqs..., trans_eqs...], t, [S, I, R, β], [γ, μ]; tspan = tspan)
ode_simple = structural_simplify(ode)

@named sde = SDESystem(
    [eqs..., trans_eqs...], noise_eqs, t, [S, I, R, β], [γ, μ, σ]; tspan = tspan
)

#%%
R₀ = 10.0
p = Dict(γ => 1 / 8, μ => 1 / (62 * 365), σ => 0.01)
u₀ = Dict(S => 999.0, I => 1.0, R => 0.0)
push!(u₀, β => calculate_beta(R₀, p[γ], p[μ], ones(1, 1), [sum(values(u₀))]))

#%%
sde_prob = SDEProblem(sde, u₀, tspan, p; check_length = false)
sde_sol = solve(sde_prob, SOSRI(); saveat = δt)
sde_sol_df, sde_sol_beta = create_sir_betas_df(sde_sol)
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
ensemble_prob = EnsembleProblem(sde_prob)
ensemble_sol = solve(ensemble_prob, EnsembleThreads(); trajectories = nsims, saveat = δt)

#%%
all_sims_array = fill(NaN, 5, tlength, nsims)
create_sir_all_sims_array!(ensemble_sol, nsims; β = true)

all_sims_array[5, :, :] = exp.(all_sims_array[5, :, :])

#%%
quantiles = [0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975]
sim_quantiles = zeros(Float64, length(quantiles), tlength, 5)
create_sir_all_sim_quantiles!(; quantiles = quantiles)

#%%
create_sir_quantiles_plot(; lower = 0.025, upper = 0.975, quantiles = quantiles)

#%%
create_beta_quantiles_plot(; lower = 0.025, upper = 0.975, quantiles = quantiles)