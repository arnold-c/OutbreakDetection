"""
This is a simulation of an SIR model with a time-dependent infectivity rate.
It uses the JumpProcesses.jl package and the Coevolve method.
It uses VariableRateJumps to simulate the infectivity rate of an individual being a function of time of the year i.e. the infectivity rate peaks once a year.
All jumps are manually defined using JumpProcesses, not using the ModelingToolkit or Catalyst interfaces.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using JumpProcesses, Statistics, DataFrames, DataFramesMeta, LinearAlgebra
using CairoMakie, AlgebraOfGraphics, DifferentialEquations, ModelingToolkit

CairoMakie.activate!()
set_aog_theme!()

#%%
# Revise will keep track of file changes and reload functions as necessary
includet(srcdir("Julia/transmission-functions.jl"))
includet(srcdir("Julia/plotting-functions.jl"))
includet(srcdir("Julia/cleaning-functions.jl"))

#%%
u₀ = [999, 1, 0]
δt = 0.5
tlower = 0.0
tmax = 365.0 * 100
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

dur_inf = 8
R₀ = 10.0
γ = 1 / dur_inf
μ = 1 / (62 * 365)
β = calculate_beta(R₀, γ, μ, 1, sum(u₀))
# Adjust the scale of the seasonal variation in infectivity i.e. A_scale scales the amplitude of cosine function to 1/A_scale. The maximum infectivity is β
A_scale = 2.0
ε = (1.06 * μ * (R₀ - 1)) / sqrt(sum(u₀)) # Commuter imports - see p210 Keeling & Rohani
p = (β, γ, μ, A_scale, ε)

#%%
birth_rate(u, p, t) = p[3] * (u[1] + u[2] + u[3])  # μ*N
birth_affect!(integrator) = integrator.u[1] += 1  # S -> S + 1
birth_jump = ConstantRateJump(birth_rate, birth_affect!)

recov_rate(u, p, t) = p[2] * u[2]         # γ*I
function recov_affect!(integrator)
    integrator.u[2] -= 1
    integrator.u[3] += 1
    return nothing
end
recov_jump = ConstantRateJump(recov_rate, recov_affect!)

S_death_rate(u, p, t) = p[3] * u[1]  # μ*S
S_death_affect!(integrator) = integrator.u[1] -= 1  # S -> S + 1
S_death_jump = ConstantRateJump(S_death_rate, S_death_affect!)

I_death_rate(u, p, t) = p[3] * u[2]  # μ*I
I_death_affect!(integrator) = integrator.u[2] -= 1  # S -> S + 1
I_death_jump = ConstantRateJump(I_death_rate, I_death_affect!)

R_death_rate(u, p, t) = p[3] * u[3]  # μ*R
R_death_affect!(integrator) = integrator.u[3] -= 1  # S -> S + 1
R_death_jump = ConstantRateJump(R_death_rate, R_death_affect!)

import_rate(u, p, t) = p[5] * (u[1] + u[2] + u[3]) / R₀   # ε*N/R₀
function import_affect!(integrator)
    integrator.u[1] -= 1    # S -> S - 1
    integrator.u[2] += 1    # I -> I + 1
    return nothing
end
import_jump = ConstantRateJump(import_rate, import_affect!)

# export_rate(u, p, t) = p[5] * (u[1] + u[2] + u[3]) / R₀ # ε
# export_affect!(integrator) = integrator.u[3] -= 1  # R -> R - 1
# export_jump = ConstantRateJump(export_rate, export_affect!)

# Place infection at the end as it is a VariableRateJump, which is ordered after ConstantRateJumps in the dependency graph.
# Amplitude = 0.5 * (cos(2pi * t / 365) + 1)) * 1/scale + (scale - 1)/scale
function infec_rate(u, p, t)
    # β*S*I * 0.5*cos(2πt/365) + 0.5
    Amplitude = 0.5 * (cos(2π * t / 365) + 1) * 1 / p[4] + (p[4] - 1) / p[4]
    return p[1] * u[1] * u[2] * Amplitude
end

# Lower bound on infection rate
lrate(u, p, t) = 0.0
# Upper bound on infection rate
urate(u, p, t) = p[1] * u[1] * u[2]  # β*S*I
# Interval that the infection rate must be within the bounds for
rateinterval(u, p, t) = 1 / (2 * urate(u, p, t))

function infec_affect!(integrator)
    integrator.u[1] -= 1     # S -> S - 1
    integrator.u[2] += 1     # I -> I + 1
    return nothing
end
infec_jump = VariableRateJump(infec_rate, infec_affect!; lrate, urate, rateinterval)

jumps = JumpSet(
    birth_jump, recov_jump, S_death_jump, I_death_jump, R_death_jump,
    import_jump,
    infec_jump,
)

# ConstantRateJumps are ordered before VariableRateJumps in the dependency graph, otherwise within the same ordering presented to the JumpProblem, i.e., birth, recovery ...
# Each row of the dependency graph is a list of rates that must be recalculated when the row's jump is executed, as it affects the underlying states each of the rates depends on.
dep_graph = [
    [1, 3, 7],          # Birth, S death, infection
    [2, 1, 4, 5, 7],    # Recovery, birth, I death, R death, infection
    [3, 1, 7],          # S death, birth, infection
    [4, 1, 2, 7],       # I death, birth, recovery, infection
    [5, 1, 2, 7],       # R death, birth, recovery, infection
    [6, 1, 2],          # Import, infection
    [7, 1, 2, 3, 4],    # Infection, birth, recovery, S death, I death
]

#%%
season_infec_prob = DiscreteProblem(u₀, tspan, p)
# Bounded VariableJumpRate problems require the Coevolve() algorithm
season_infec_jump_prob = JumpProblem(season_infec_prob, Coevolve(), jumps; dep_graph)
season_infec_sol = solve(season_infec_jump_prob, SSAStepper())

season_infec_sol_df = create_sir_df(season_infec_sol, [:S, :I, :R])

#%%
sircolors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

draw_sir_plot(season_infec_sol_df; colors = sircolors, annual = true)

#%%
# Visualize seasonal changes in the infection rate
t = 0:0.01:720
y = Vector{Float64}(undef, length(t))
y_scale = similar(y)
@. y = 0.5 * (cos(2pi * t / 365) + 1)
@. y_scale = (0.5 * (cos(2pi * t / 365) + 1)) * 1 / A_scale + (A_scale - 1) / A_scale
fig = lines(t, y; color = "black", linewidth = 2)
lines!(t, y_scale; color = "green", linewidth = 2)
vlines!(0:365:720; color = (:red, 0.5), linestyle = :dash, linewidth = 2)
hlines!([0.0, 0.5, 1.0])
fig

#%%
nsims = 1000

sir_array = zeros(4, tlength)
all_sims_array = fill(NaN, 4, tlength, nsims)

season_infec_jump_prob = JumpProblem(
    season_infec_prob, Coevolve(), jumps; dep_graph = dep_graph,
    save_positions = (false, false),
)

create_sir_all_sims_array!(;
    nsims = nsims, prob = season_infec_jump_prob, alg = SSAStepper(), δt = δt
)

#%%
quantiles = [0.025, 0.05, 0.1, 0.2, 0.25, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9, 0.95, 0.975]
sim_quantiles = zeros(Float64, length(quantiles), tlength, 4)

#%%
create_sir_all_sim_quantiles!(; quantiles = quantiles)

#%%
create_sir_quantiles_plot(; lower = 0.1, upper = 0.9, quantiles = quantiles)