"""
This is a stochastic differential equation SIR model using the ModelingToolkit interface.
The final section of code calculates the R₀ value, as well as the β parameter, for the model, using the Next Generation Matrix method.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using DifferentialEquations, ModelingToolkit, Statistics
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
eqs = [D(S) ~ -β * S * I + μ * (I + R), D(I) ~ β * S * I - γ * I - μ * I,
    D(R) ~ γ * I - μ * R]

noise_eqs = [
    0.5 + 0.01 * S,
    0.5 + 0.01 * I,
    0.5 + 0.01 * R,
]

#%%
@named ode = ODESystem(eqs, t, [S, I, R], [β, γ, μ]; tspan = tspan)
ode_simple = structural_simplify(ode)

@named sde = SDESystem(ode, noise_eqs)

#%%
R₀ = 10.0
p = Dict(γ => 1 / 8, μ => 1 / (62 * 365))
u₀ = Dict(S => 999.0, I => 1.0, R => 0.0)
push!(
    p,
    β => calculate_beta(ode_simple, 1, 1, R₀, p, [1.0], [sum(values(u₀))]),
)

#%%
sde_prob = SDEProblem(sde, u₀, tspan, p)
sde_sol = solve(sde_prob, SOSRI(); saveat = δt)
sde_sol_df = create_sir_df(sde_sol)

#%%
colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(sde_sol_df; colors = colors)

#%%
nsims = 1000

sir_array = zeros(5, tlength)
all_sims_array = fill(NaN, 5, tlength, nsims)

quantiles = [0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975]
sim_quantiles = zeros(Float64, length(quantiles), tlength, 4)


#%%
ensemble_prob = EnsembleProblem(sde_prob)
ensemble_sol = solve(ensemble_prob, EnsembleThreads(); trajectories = nsims, saveat = δt)

ensemble_summ = EnsembleSummary(ensemble_sol, tlower:1:tmax; quantiles = [0.025, 0.975])

create_sir_plot(create_sir_df(ensemble_summ.qhigh))

function create_sir_ensemble_df(ensemble::EnsembleSummary, quantiles::Vector{Float64})
    med_df = transform!(create_sir_df(ensemble.med), y -> (y = "median"))
    low_df = transform!(create_sir_df(ensemble.qlow), y -> (y = "lower"))
    high_df = transform!(create_sir_df(ensemble.qhigh), y -> (y = "upper"))

    return vcat(med_df, low_df, high_df)
end

sim_quantiles = zeros(4, length(ensemble_summ.med), length(quantiles))
perm_sim_quantiles = permutedims(sim_quantiles, (3, 2, 1))


function create_sir_ensemble_array!(
    ensemble::EnsembleSummary; lower::T, upper::T, quantiles::Vector{Float64}
) where {(T <: AbstractFloat)}
    med_index = findfirst(isequal(0.5), quantiles)
    lower_index = findfirst(isequal(lower), quantiles)
    upper_index = findfirst(isequal(upper), quantiles)

    sim_quantiles[1:3, :, med_index] .= Array(ensemble.med)
    sim_quantiles[1:3, :, lower_index] = Array(ensemble.qlow)
    sim_quantiles[1:3, :, upper_index] = Array(ensemble.qhigh)

    sim_quantiles[4, :, :] = sum(sim_quantiles[1:3, :, :]; dims=1)
    permutedims!(perm_sim_quantiles, sim_quantiles, (3, 2, 1))
    
    return perm_sim_quantiles
end

create_sir_ensemble_array!(ensemble_summ; lower = 0.025, upper = 0.975, quantiles = quantiles)

create_sir_quantiles_plot(perm_sim_quantiles; lower = 0.025, upper = 0.975, quantiles = quantiles, δt = 1.0)
