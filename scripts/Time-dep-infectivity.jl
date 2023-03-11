"""
This is a simulation of an SIR model with a time-dependent infectivity rate.
It uses the JumpProcesses.jl package and the Coevolve method.
It uses VariableRateJumps to simulate the infectivity rate of an individual being a function of time since infection i.e. the infectivity rate decays from a peak infectivity of β₁ + α exponentially until the individual recovers.
All jumps are manually defined using JumpProcesses, not using the ModelingToolkit or Catalyst interfaces.
"""

using DrWatson
@quickactivate "OutbreakDetection"

using JumpProcesses, Statistics, DataFrames, DataFramesMeta, LinearAlgebra
using CairoMakie, AlgebraOfGraphics

CairoMakie.activate!()
set_aog_theme!()

includet(srcdir("Julia/transmission-functions.jl")) # Revise will keep track of file changes and reload functions as necessary
includet(srcdir("Julia/plotting-functions.jl"))
includet(srcdir("Julia/cleaning-functions.jl"))

u₀ = [999, 1, 0]
tspan = (0.0, 250)

dur_inf = 8
R₀ = 2.0
γ = 1/dur_inf
μ = 1/(62 * 365)
α =  calculate_beta(R₀, γ, μ, 1, sum(u₀))
β₁ = α / 1000
ν = 0.05 # Decay rate
p = (β₁, γ, μ, α, ν)

# Holds the timestamp of active infections, therefore is the same length as the duration of infection at the start of the simulation. Each infection adds a timestamp to the end of the array and each recovery removes a timestamp from a random position of the array.
H = zeros(Float64, dur_inf)

birth_rate(u, p, t) = p[3] * (u[1] + u[2] + u[3])  # μ*N
birth_affect!(integrator) = integrator.u[1] += 1  # S -> S + 1
birth_jump = ConstantRateJump(birth_rate, birth_affect!)

recov_rate(u, p, t) = p[2] * u[2]         # γ*I
function recov_affect!(integrator)
    integrator.u[2] -= 1
    integrator.u[3] += 1
    length(H) > 0 && deleteat!(H, rand(1:length(H)))
    nothing
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

# Place infection at the end as it is a VariableRateJump, which is ordered after ConstantRateJumps in the dependency graph.
function infec_rate(u, p, t)
    # β₁*S*I + α*S*∑ₜ exp(-ν(t - tᵢ))
    p[1] * u[1] * u[2] + p[4] * u[1] * sum(exp(-p[5] * (t - _t)) for _t in H)
end

# Lower bound on infection rate
lrate(u, p, t) = p[1] * u[1] * u[2]  # β₁*S*I
# Upper bound on infection rate
urate = infec_rate
# Interval that the infection rate must be within the bounds for
rateinterval(u, p, t) = 1 / (2 * urate(u, p, t))

function infec_affect!(integrator)
    integrator.u[1] -= 1     # S -> S - 1
    integrator.u[2] += 1     # I -> I + 1
    push!(H, integrator.t)
    nothing
end
infec_jump = VariableRateJump(infec_rate, infec_affect!; lrate, urate, rateinterval)

jumps = [birth_jump, recov_jump, S_death_jump, I_death_jump, R_death_jump, infec_jump]

# ConstantRateJumps are ordered before VariableRateJumps in the dependency graph, otherwise within the same ordering presented to the JumpProblem, i.e., birth, recovery ...
# Each row of the dependency graph is a list of rates that must be recalculated when the row's jump is executed, as it affects the underlying states each of the rates depends on.
dep_graph = [
    [1, 3, 6],          # Birth, S death, infection
    [2, 1, 4, 5, 6],    # Recovery, birth, I death, R death, infection
    [3, 1, 6],          # S death, birth, infection
    [4, 1, 2, 6],       # I death, birth, recovery, infection
    [5, 1, 2, 6],       # R death, birth, recovery, infection
    [6, 1, 2, 3, 4]     # Infection, birth, recovery, S death, I death
    ]

t_dep_infec_prob = DiscreteProblem(u₀, tspan, p)
# Bounded VariableJumpRate problems require the Coevolve() algorithm
t_dep_infec_jump_prob = JumpProblem(
    t_dep_infec_prob, Coevolve(), jumps...; dep_graph
    )
t_dep_infec_sol = solve(t_dep_infec_jump_prob, SSAStepper())

t_dep_infec_sol_df = create_sir_df(t_dep_infec_sol)

colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(t_dep_infec_sol_df; colors = colors)