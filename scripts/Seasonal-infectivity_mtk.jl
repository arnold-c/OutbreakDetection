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
tlength = length(tlower:δt:tmax);

#%%
@parameters β γ μ A_scale t
@variables S(t) I(t) R(t)

#%%
u₀ = [999.0, 1.0, 0.0]
R₀ = 10.0
gamma = 1 / 8
mu = 1 / (62 * 365)
amp_scale = 2.0
beta = calculate_beta(R₀, gamma, mu, ones(1, 1), [sum(u₀)])
p = (beta, gamma, mu, amp_scale)

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
# A = 0.5 * (cos(2π * t / 365) + 1) * 1 / A_scale + (A_scale - 1) / A_scale

# Lower bound on infection rate
lrate = 0.0
# Upper bound on infection rate
urate = β * S * I
# Interval that the infection rate must be within the bounds for
rateinterval = 1 / (2 * urate)
# * (0.5 * (cos(2π * t / 365) + 1) * 1 / A_scale + (A_scale - 1) / A_scale)
infec_rate = β * S * I * (0.5 * (cos(2π * t / 365) + 1) * 1 / A_scale + (A_scale - 1) / A_scale)
infec_affect! = [S ~ S - 1, I ~ I + 1]

#%%
recov_jump = ConstantRateJump(recov_rate, recov_affect!)
birth_jump = ConstantRateJump(birth_rate, birth_affect!)
S_death_jump = ConstantRateJump(S_death_rate, S_death_affect!)
I_death_jump = ConstantRateJump(I_death_rate, I_death_affect!)
R_death_jump = ConstantRateJump(R_death_rate, R_death_affect!)
infec_jump = VariableRateJump(infec_rate, infec_affect!; lrate = lrate, urate = urate, rateinterval = rateinterval)

jumps = [recov_jump, birth_jump, S_death_jump, I_death_jump, R_death_jump, infec_jump]

#%%
@named js = JumpSystem(jumps, t, [S, I, R], [β, γ, μ, A_scale])
ode_prob = ODEProblem((du, u, p, t) -> du .= 0, u₀, tspan, p)
js_prob = JumpProblem(js, ode_prob, Direct())
js_sol = solve(js_prob, Tsit5())
js_sol_df = create_sir_df(js_sol)

#%%
colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(js_sol_df; colors = colors)

#%%
# Visualize seasonal changes in the infection rate
t = 0:0.01:720
y = Vector{Float64}(undef, length(t))
y_scale = similar(y)
@. y = 0.5 * (cos(2pi * t / 365) + 1)
@. y_scale = (0.5 * (cos(2pi * t / 365) + 1)) * 1 / amp_scale + (amp_scale - 1) / amp_scale
fig = lines(t, y; color = "black", linewidth = 2)
lines!(t, y_scale; color = "green", linewidth = 2)
vlines!(0:365:720; color = (:red, 0.5), linestyle = :dash, linewidth = 2)
hlines!([0.0, 0.5, 1.0])
fig

#%%
nsims = 1000

sir_array = zeros(4, tlength)
all_sims_array = fill(NaN, 4, tlength, nsims)

season_infec_jump_prob =  JumpProblem(js, ode_prob, Direct(); save_positions = (false, false))

create_sir_all_sims_array!(;
    nsims = nsims, prob = season_infec_jump_prob, alg = Tsit5(), δt = δt
)

solve(season_infec_jump_prob, Tsit5(); saveat = δt)

#%%
quantiles = [0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 0.8, 0.9, 0.95, 0.975]
sim_quantiles = zeros(Float64, length(quantiles), tlength, 4)

#%%
create_sir_all_sim_quantiles!(; quantiles = quantiles)

#%%
create_sir_quantiles_plot(; lower = 0.1, upper = 0.9, quantiles = quantiles)

# #%%
# @parameters β γ t
# @variables S(t) I(t) R(t)

# rate₁ = β * S * I
# affect₁ = [S ~ S - 1, I ~ I + 1]
# rate₂ = γ * I * exp(-t / 20)
# affect₂ = [I ~ I - 1, R ~ R + 1]

# j₁ = ConstantRateJump(rate₁, affect₁)
# j₂ = VariableRateJump(rate₂, affect₂)

# @named js = JumpSystem([j₁, j₂], t, [S, I, R], [β, γ])

# u₀ = [999.0, 1.0, 0.0];
# p = (0.1 / 1000, 0.01);
# tspan = (0.0, 250.0);

# oprob = ODEProblem((du, u, p, t) -> du .= 0, u₀, tspan, p)
# jprob = JumpProblem(js, oprob, Direct())
# sol = solve(jprob, Tsit5())