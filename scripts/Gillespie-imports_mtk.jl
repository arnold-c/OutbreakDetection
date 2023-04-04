"""
This is a Gillespie simulation of an SIR model that has imported infections.
It uses the JumpProcesses.jl package and the Direct Gillespie method.
All jumps are constant rate jumps and defined using the ModelingToolkit interface.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using DifferentialEquations, ModelingToolkit, Statistics
using DataFrames, DataFramesMeta, LinearAlgebra
using CairoMakie, AlgebraOfGraphics

CairoMakie.activate!()
set_aog_theme!()

#%%
includet(srcdir("Julia/transmission-functions.jl")) # Revise will keep track of file changes and reload functions as necessary
includet(srcdir("Julia/plotting-functions.jl"))
includet(srcdir("Julia/cleaning-functions.jl"))

#%%
@parameters β γ μ σ
@variables t S(t) I(t) R(t)

u₀ = Dict(S => 999, I => 1, R => 0)
R₀ = 10.0

#%%
p = Dict(γ => 1 / 8, μ => 1 / (62 * 365), σ => 1 / (365 * 10))
push!(p, β => calculate_beta(R₀, p[μ], p[γ], 1, sum(values(u₀))))

#%%
δt = 0.005
tlower = 0.0
tmax = 365.0 * 100
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

#%%
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
import_rate = σ * (S + I + R)
import_affect! = [I ~ I + 1]
export_rate = import_rate
export_affect! = [R ~ R - 1]

#%%
infec_jump = ConstantRateJump(infec_rate, infec_affect!)
recov_jump = ConstantRateJump(recov_rate, recov_affect!)
birth_jump = ConstantRateJump(birth_rate, birth_affect!)
S_death_jump = ConstantRateJump(S_death_rate, S_death_affect!)
I_death_jump = ConstantRateJump(I_death_rate, I_death_affect!)
R_death_jump = ConstantRateJump(R_death_rate, R_death_affect!)
import_jump = ConstantRateJump(import_rate, import_affect!)
export_jump = ConstantRateJump(export_rate, export_affect!)

jumps = [infec_jump, recov_jump, birth_jump, S_death_jump, I_death_jump, R_death_jump, import_jump, export_jump]

#%%
@named js = JumpSystem(jumps, t, [S, I, R], [β, γ, μ, σ])

gill_dis_prob = DiscreteProblem(js, u₀, tspan, p)
gill_jump_prob = JumpProblem(js, gill_dis_prob, Direct())
gill_sol = solve(gill_jump_prob, SSAStepper())

gill_sol_df = create_sir_df(gill_sol)

#%%
colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(gill_sol_df; colors = colors)

#%%
nsims = 1000

gill_jump_prob = JumpProblem(js, gill_dis_prob, Direct(); save_positions = (false, false))

#%%
sir_array = zeros(4, tlength)
all_sims_array = fill(NaN, 4, tlength, nsims)

create_sir_all_sims_array!(;
nsims = nsims, prob = gill_jump_prob, alg = SSAStepper(), δt = δt
)

#%%
quantiles = [0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 0.8, 0.9, 0.95, 0.975]
sim_quantiles = zeros(Float64, length(quantiles), tlength, 4)

create_sir_all_sim_quantiles!(; quantiles = quantiles)

#%%
create_sir_quantiles_plot(; lower = 0.2, upper = 0.8, quantiles = quantiles)