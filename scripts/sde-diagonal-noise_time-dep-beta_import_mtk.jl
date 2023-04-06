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
tmax = 365.0 * 100
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

#%%
@parameters γ μ σ ε
@variables t S(t) I(t) R(t) β(t)
D = Differential(t)

#%%
nic = 1 # number of infective compartments
nac = 1 # number of age classes

#%%
eqs = [
    D(S) ~ -exp(β) * S * I + μ * (I + R),
    D(I) ~ exp(β) * S * I - γ * I - μ * I + ε,
    D(R) ~ γ * I - μ * R - ε,
]

trans_eqs = [D(β) ~ 0.0]

noise_eqs = [
    0.0,
    0.0,
    0.0,
    0.0,
]

#%%
@named ode = ODESystem([eqs..., trans_eqs...], t, [S, I, R, β], [γ, μ, ε]; tspan = tspan)
ode_simple = structural_simplify(ode)

@named sde = SDESystem(
    [eqs..., trans_eqs...], noise_eqs, t, [S, I, R, β], [γ, μ, σ, ε]; tspan = tspan
)

#%%
R₀ = 10.0
p = Dict(γ => 1 / 8, μ => 1 / (62 * 365), σ => 0.01, ε => 0.0)
u₀ = Dict(S => 999.0, I => 1.0, R => 0.0)

logbeta = log(calculate_beta(R₀, p[γ], p[μ], ones(1, 1), [sum(values(u₀))]))
push!(u₀, β => logbeta)

#%%
sde_prob = SDEProblem(sde, u₀, tspan, p; check_length = false)
sde_sol = solve(sde_prob, SOSRA(); saveat = δt)
sde_sol_df, sde_sol_beta = create_sir_beta_dfs(sde_sol, [:S, :I, :R])
sde_sol_beta.β = exp.(sde_sol_beta.β)

#%%
beta_plot = create_beta_plot(sde_sol_beta)
draw_beta_plot(beta_plot)

#%%
sircolors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

sir_plot = create_sir_plot(sde_sol_df)
draw_sir_plot(sir_plot; colors = sircolors)

#%%
draw_combined_sir_beta_plot(sir_plot, beta_plot)

#%%
ode_prob = ODEProblem(ode_simple, u₀, tspan, p)
ode_sol = solve(ode_prob, Tsit5(); saveat = δt)
ode_sol_df, ode_sol_beta = create_sir_beta_dfs(ode_sol, [:S, :I, :R])
ode_sol_beta.β = exp.(ode_sol_beta.β)

draw_beta_plot(ode_sol_beta)
draw_sir_plot(ode_sol_df; colors = sircolors)

@chain ode_sol_df begin
    subset(_, :Number => num -> num .>= 0.0 .&& num .<= 1000.0)
    draw_sir_plot(_; colors = sircolors)
end

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