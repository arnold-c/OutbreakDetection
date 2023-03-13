"""
This is a basic Gillespie simulation of an SIR model.
It uses the JumpProcesses.jl package and the Direct Gillespie method.
All jumps are constant rate jumps and are manually defined using JumpProcesses, not using the ModelingToolkit or Catalyst interfaces.
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
δt = 0.005
tlower = 0.0
tmax = 250.0
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

R₀ = 10.0
γ = 1/8
μ = 1/(62 * 365)

β = calculate_beta(R₀, γ, μ, 1, sum(u₀))
p = (β, γ, μ)

birth_rate(u, p, t) = p[3] * (u[1] + u[2] + u[3])  # μ*N
birth_affect!(integrator) = integrator.u[1] += 1  # S -> S + 1
birth_jump = ConstantRateJump(birth_rate, birth_affect!)

infec_rate(u, p, t) = p[1] * u[1] * u[2]  # β*S*I
function infec_affect!(integrator)
    integrator.u[1] -= 1         # S -> S - 1
    integrator.u[2] += 1         # I -> I + 1
    nothing
end
infec_jump = ConstantRateJump(infec_rate, infec_affect!)

recov_rate(u, p, t) = p[2] * u[2]         # ν*I
function recov_affect!(integrator)
    integrator.u[2] -= 1        # I -> I - 1
    integrator.u[3] += 1        # R -> R + 1
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

jumps = [birth_jump, infec_jump, recov_jump, S_death_jump, I_death_jump, R_death_jump]

gill_prob = DiscreteProblem(u₀, tspan, p)

gill_jump_prob = JumpProblem(gill_prob, Direct(), jumps...)
gill_sol = solve(gill_jump_prob, SSAStepper())

gill_sol_df = create_sir_df(gill_sol)

colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(gill_sol_df, colors = colors)

nsims = 1000

sir_array = zeros(5, tlength)

all_sims_array = fill(NaN, 5, tlength, nsims)

quantiles = [0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975]

sim_quantiles = zeros(Float64, length(quantiles), tlength, 4)

gill_jump_prob = JumpProblem(
    gill_prob, Direct(), jumps...; save_positions = (false, false)
    )


create_sir_all_sims_array!(;
    nsims = nsims, prob = gill_jump_prob, alg = SSAStepper(), δt = δt
    )

create_sir_all_sim_quantiles!(quantiles = quantiles)

create_sir_quantiles_plot!(lower = 0.1, upper = 0.9, quantiles = quantiles)