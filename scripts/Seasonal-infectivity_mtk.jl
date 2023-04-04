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
δt = 0.005
tlower = 0.0
tmax = 250.0
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

#%%
@parameters β γ μ A_scale
@variables t S(t) I(t) R(t)
D = Differential(t)

#%%
nic = 1 # number of infective compartments
nac = 1 # number of age classes

#%%
# Create ODE system for β calculation using NGM
ode_eqs = [D(S) ~ -β * S * I + μ * (I + R), D(I) ~ β * S * I - γ * I - μ * I,
    D(R) ~ γ * I - μ * R]

@named de = ODESystem(ode_eqs, t, [S, I, R], [β, γ, μ]; tspan = tspan)
de_simple = structural_simplify(de)

#%%
R₀ = 10.0
p = Dict(γ => 1 / 8, μ => 1 / (62 * 365), A_scale => 10.0)
u₀ = Dict(S => 999.0, I => 1.0, R => 0.0)
push!(p, β => calculate_beta(de_simple, 1, 1, R₀, p, [1.0], [sum(values(u₀))]))

#%%
recov_rate = γ * I
recov_affect! = [I ~ I - 1, R ~ R + 1]
birth_rate = μ * (S + I + R)
birth_affect! = [S ~ S + 1]
S_death_rate = μ * S
S_death_affect! = [S ~ S - 1]
I_death_rate = μ * I
I_death_affect! = [I ~ I - 1]
R_death_rate = μ * R
R_death_affect! = [R ~ R - 1]

# Place infection at the end as it is a VariableRateJump, which is ordered after ConstantRateJumps in the dependency graph.
# A = 0.5 * (cos(2pi * t / 365) + 1)) * 1/scale + (scale - 1)/scale
A = 0.5 * (cos(2π * t / 365) + 1) * 1 / A_scale + (A_scale - 1) / A_scale

# Lower bound on infection rate
lrate = 0.0
# Upper bound on infection rate
urate = β * S * I
# Interval that the infection rate must be within the bounds for
rateinterval = 1 / (2 * urate)

infec_rate = β * S * I * A
infec_affect! = [S ~ S - 1, I ~ I + 1]

#%%
recov_jump = ConstantRateJump(recov_rate, recov_affect!)
birth_jump = ConstantRateJump(birth_rate, birth_affect!)
S_death_jump = ConstantRateJump(S_death_rate, S_death_affect!)
I_death_jump = ConstantRateJump(I_death_rate, I_death_affect!)
R_death_jump = ConstantRateJump(R_death_rate, R_death_affect!)
infec_jump = VariableRateJump(infec_rate, infec_affect!; lrate, urate, rateinterval)

jumps = [recov_jump, birth_jump, S_death_jump, I_death_jump, R_death_jump, infec_jump]

#%%
# ConstantRateJumps are ordered before VariableRateJumps in the dependency graph, otherwise within the same ordering presented to the JumpProblem, i.e., birth, recovery ...
# Each row of the dependency graph is a list of rates that must be recalculated when the row's jump is executed, as it affects the underlying states each of the rates depends on.
dep_graph = [
    [1, 3, 6],          # Birth, S death, infection
    [2, 1, 4, 5, 6],    # Recovery, birth, I death, R death, infection
    [3, 1, 6],          # S death, birth, infection
    [4, 1, 2, 6],       # I death, birth, recovery, infection
    [5, 1, 2, 6],       # R death, birth, recovery, infection
    [6, 1, 2, 3, 4],     # Infection, birth, recovery, S death, I death
]

#%%
@named js = JumpSystem(jumps, t, [S, I, R], [β, γ, μ, A_scale]; dep_graph)

#%%
de_prob = DiscreteProblem(js, u₀, tspan, p)
js_prob = JumpProblem(de_prob, Coevolve())
js_sol = solve(js_prob, SSAStepper())
js_sol_df = create_sir_df(js_sol)

#%%
season_infec_prob = DiscreteProblem([values(u₀)...], tspan, [values(p)...])
# Bounded VariableJumpRate problems require the Coevolve() algorithm
season_infec_jump_prob = JumpProblem(season_infec_prob, Coevolve(), jumps...; dep_graph)
season_infec_sol = solve(season_infec_jump_prob, SSAStepper())

season_infec_sol_df = create_sir_df(season_infec_sol)

#%%
colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(js_sol_df; colors = colors)

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
    season_infec_prob, Coevolve(), jumps...; dep_graph, save_positions = (false, false)
)

create_sir_all_sims_array!(;
    nsims = nsims, prob = season_infec_jump_prob, alg = SSAStepper(), δt = δt
)

#%%
quantiles = [0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 0.8, 0.9, 0.95, 0.975]
sim_quantiles = zeros(Float64, length(quantiles), tlength, 4)

#%%
create_sir_all_sim_quantiles!(; quantiles = quantiles)

#%%
create_sir_quantiles_plot(; lower = 0.1, upper = 0.9, quantiles = quantiles)