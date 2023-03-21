"""
This is a Gillespie simulation of an SIR model that has imported infections.
It uses the JumpProcesses.jl package and the Direct Gillespie method.
All jumps are constant rate jumps and defined using the ModelingToolkit interface.
"""

using DrWatson
@quickactivate "OutbreakDetection"

using DifferentialEquations, ModelingToolkit, Statistics
using DataFrames, DataFramesMeta, LinearAlgebra, Symbolics
using CairoMakie, AlgebraOfGraphics

CairoMakie.activate!()
set_aog_theme!()

includet(srcdir("Julia/transmission-functions.jl")) # Revise will keep track of file changes and reload functions as necessary
includet(srcdir("Julia/plotting-functions.jl"))
includet(srcdir("Julia/cleaning-functions.jl"))

u₀ = Dict(S => 999, I => 1, R => 0)
R₀ = 10.0

p = Dict(γ => 1 / 8, μ => 1 / (62 * 365))
push!(p, β => calculate_beta(R₀, p[μ], p[γ], 1, sum(values(u₀))))

δt = 0.005
tlower = 0.0
tmax = 250.0
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

@parameters β γ μ
@variables t S(t) I(t) R(t)
D = Differential(t)

eqs = [
    D(S) ~ -β * S * I + μ * (I + R), D(I) ~ β * S * I - γ * I - μ * I, D(R) ~ γ * I - μ * R
]

@named de = ODESystem(eqs, t, [S, I, R], [β, γ, μ]; tspan=tspan)
de_simple = structural_simplify(de)

ode_prob = ODEProblem(de_simple, u₀, tspan, p)
ode_sol = solve(ode_prob, Tsit5(); saveat=δt)

nic = 1 # number of infective compartments
nac = 1 # number of age classes

Jac = calculate_jacobian(de_simple)[
    (nac + 1):(nac + nic * nac), (nac + 1):(nac + nic * nac)
]

F = substitute(Jac, Dict(γ => 0.0, μ => 0.0))
V = substitute(Jac, Dict(β => 0.0))
FV⁻¹ = F * -inv(V)
all_eigenvals = eigvals(eigen(substitute(FV⁻¹, Dict(S => 1000.0, p...))))
R0 = maximum(real(all_eigenvals))

eigen(
    substitute(
        calculate_jacobian(de_simple)[
            (nac + 1):(nac + nic * nac), (nac + 1):(nac + nic * nac)
        ],
        Dict(S => 0.0, I => 0.0, R => 0.0, p...),
    ),
)

gill_dis_prob = DiscreteProblem(js, u₀, tspan, p)
gill_jump_prob = JumpProblem(js, gill_dis_prob, Direct(); jac=true)
gill_sol = solve(gill_jump_prob, SSAStepper())

gill_sol_df = create_sir_df(gill_sol)

colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(gill_sol_df; colors=colors)

nsims = 1000

sir_array = zeros(5, tlength)

all_sims_array = fill(NaN, 5, tlength, nsims)

quantiles = [0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975]

sim_quantiles = zeros(Float64, length(quantiles), tlength, 4)

gill_jump_prob = JumpProblem(js, gill_dis_prob, Direct(); save_positions=(false, false))

create_sir_all_sims_array!(; nsims=nsims, prob=gill_jump_prob, alg=SSAStepper(), δt=δt)

create_sir_all_sim_quantiles!(; quantiles=quantiles)

create_sir_quantiles_plot!(; lower=0.1, upper=0.9, quantiles=quantiles)
