"""
This is a Gillespie simulation of an SIR model that has imported infections.
It uses the JumpProcesses.jl package and the Direct Gillespie method.
All jumps are constant rate jumps and defined using the ModelingToolkit interface.
"""

using DrWatson
@quickactivate "OutbreakDetection"

using DifferentialEquations, ModelingToolkit, Statistics
using DataFrames, DataFramesMeta, LinearAlgebra
using CairoMakie, AlgebraOfGraphics

CairoMakie.activate!()
set_aog_theme!()

includet(srcdir("Julia/transmission-functions.jl")) # Revise will keep track of file changes and reload functions as necessary
includet(srcdir("Julia/plotting-functions.jl"))
includet(srcdir("Julia/cleaning-functions.jl"))

u₀ = Dict(S => 999, I => 1, R => 0)
R₀ = 10.0

p = Dict(γ => 1/8, μ => 1/(62 * 365))
push!(p, β => calculate_beta(R₀, p[μ], p[γ], 1, sum(values(u₀))))

δt = 0.005
tlower = 0.0
tmax = 250.0
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

@parameters β γ μ
@variables t S(t) I(t) R(t)

infec_rate = β * S * I
infec_affect! = [S ~ S - 1, I ~ I + 1]
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

infec_jump = ConstantRateJump(infec_rate, infec_affect!)
recov_jump = ConstantRateJump(recov_rate, recov_affect!)
birth_jump = ConstantRateJump(birth_rate, birth_affect!)
S_death_jump = ConstantRateJump(S_death_rate, S_death_affect!)
I_death_jump = ConstantRateJump(I_death_rate, I_death_affect!)
R_death_jump = ConstantRateJump(R_death_rate, R_death_affect!)

jumps = [infec_jump, recov_jump, birth_jump, S_death_jump, I_death_jump, R_death_jump]

@named js = JumpSystem(jumps, t, [S, I, R], [β, γ, μ])


gill_dis_prob = DiscreteProblem(js, u₀, tspan, p)
gill_jump_prob = JumpProblem(js, gill_dis_prob, Direct())
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
    js, gill_dis_prob, Direct(); save_positions = (false, false)
    )


create_sir_all_sims_array!(;
    nsims = nsims, prob = gill_jump_prob, alg = SSAStepper(), δt = δt
    )

create_sir_all_sim_quantiles!(quantiles = quantiles)

create_sir_quantiles_plot!(lower = 0.1, upper = 0.9, quantiles = quantiles)