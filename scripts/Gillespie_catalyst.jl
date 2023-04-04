"""
This is a stochastic SIR model with imports using the ModelingToolkit interface.
The final section of code calculates the R₀ value, as well as the β parameter, for the model, using the Next Generation Matrix method.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using DifferentialEquations, Catalyst, JumpProcesses, Statistics
using DataFrames, DataFramesMeta, LinearAlgebra, Symbolics
using CairoMakie, AlgebraOfGraphics

CairoMakie.activate!()
set_aog_theme!()

#%%
includet(srcdir("Julia/transmission-functions.jl")) # Revise will keep track of file changes and reload functions as necessary
includet(srcdir("Julia/plotting-functions.jl"))
includet(srcdir("Julia/cleaning-functions.jl"))

#%%
δt = 0.005
tlower = 0.0
tmax = 250.0
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

#%%
R₀ = 10.0
γ = 1 / 8
μ = 1 / (62 * 365)
S_init = 999
I_init = 1
R_init = 0

u₀ = [:S => S_init, :I => I_init, :R => R_init]

β = calculate_beta(R₀, γ, μ, 1, S_init + I_init + R_init)
p = [:γ => γ, :μ => μ, :β => β]

#%%
sir_model = @reaction_network begin
    β, S + I --> 2I
    γ, I --> R
    μ * (S + I + R), 0 --> S
    μ, (S, I, R) --> 0
end

#%%
de_prob = DiscreteProblem(sir_model, u₀, tspan, p)
jump_prob = JumpProblem(sir_model, de_prob, Direct())

#%%
jump_sol = solve(jump_prob, SSAStepper(); saveat = δt)
jump_sol_df = create_sir_df(jump_sol)

#%%
colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(jump_sol_df; colors = colors)

#%%
nsims = 1000
gill_jump_prob = JumpProblem(sir_model, de_prob, Direct(); save_positions = (false, false))

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