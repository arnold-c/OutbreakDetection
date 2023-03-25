"""
This is a stochastic SIR model with imports using the ModelingToolkit interface.
The final section of code calculates the R₀ value, as well as the β parameter, for the model, using the Next Generation Matrix method.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using DifferentialEquations, ModelingToolkit, JumpProcesses, Statistics
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
@parameters β γ μ
@variables t S(t) I(t) R(t)
D = Differential(t)

#%%
nic = 1 # number of infective compartments
nac = 1 # number of age classes

#%%
# Create ODE system for β calculation using NGM
eqs = [D(S) ~ -β * S * I + μ * (I + R), D(I) ~ β * S * I - γ * I - μ * I,
    D(R) ~ γ * I - μ * R]

@named de = ODESystem(eqs, t, [S, I, R], [β, γ, μ]; tspan = tspan)
de_simple = structural_simplify(de)

#%%
R₀ = 10.0
p = Dict(γ => 1 / 8, μ => 1 / (62 * 365))
u₀ = Dict(S => 999.0, I => 1.0, R => 0.0)
push!(p, β => calculate_beta(de_simple, 1, 1, R₀, p, [1.0], [sum(values(u₀))]))

#%%
birth_rate = μ * (S + I + R)
birth_affect = [S ~ S + 1]
birth_jump = ConstantRateJump(birth_rate, birth_affect)

infec_rate = β * S * I
infec_affect = [S ~ S - 1, I ~ I + 1]
infec_jump = ConstantRateJump(infec_rate, infec_affect)

recov_rate = γ * I
recov_affect = [I ~ I - 1, R ~ R + 1]
recov_jump = ConstantRateJump(recov_rate, recov_affect)

S_death_rate = μ * S
S_death_affect = [S ~ S - 1]
S_death_jump = ConstantRateJump(S_death_rate, S_death_affect)

I_death_rate = μ * I
I_death_affect = [I ~ I - 1]
I_death_jump = ConstantRateJump(I_death_rate, I_death_affect)

R_death_rate = μ * R
R_death_affect = [R ~ R - 1]
R_death_jump = ConstantRateJump(R_death_rate, R_death_affect)

jumps = [birth_jump, infec_jump, recov_jump, S_death_jump, I_death_jump, R_death_jump]

@named js = JumpSystem(jumps, t, [S, I, R], [β, γ, μ]; tspan = tspan)

#%%
u₀ = [S => 999, I => 1, R => 0]
p = collect([p...])

p

de_prob = DiscreteProblem(js, u₀, tspan, p)
js_prob = JumpProblem(js, de_prob, Direct())
js_sol = solve(js_prob, SSAStepper())
js_sol_df = create_sir_df(js_sol)

colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(js_sol_df; colors = colors)

#%%
nsims = 1000

sir_array = zeros(5, tlength)
all_sims_array = fill(NaN, 5, tlength, nsims)

quantiles = [0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975]
sim_quantiles = zeros(Float64, length(quantiles), tlength, 4)

#%%
gill_jump_prob = JumpProblem(js, de_prob, Direct(); save_positions = (false, false))

create_sir_all_sims_array!(;
    nsims = nsims, prob = gill_jump_prob, alg = SSAStepper(), δt = δt
)

create_sir_all_sim_quantiles!(; quantiles = quantiles)

create_sir_quantiles_plot(; lower = 0.1, upper = 0.9, quantiles = quantiles)