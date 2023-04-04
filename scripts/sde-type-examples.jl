using DifferentialEquations
using DataFrames, Plots

# Setup
f(du, u, p, t) = (du .= u)
g(du, u, p, t) = (du .= u)
u0 = rand(4, 2)

# System of SDEs with Diagonal Noise
diag_prob = SDEProblem(f, g, u0, (0.0, 1.0))
diag_sol = solve(diag_prob)
Plots.plot(diag_sol)

# System of SDEs with Scalar Noise
W = WienerProcess(0.0, 0.0, 0.0)
scalar_prob = SDEProblem(f, g, u0, (0.0, 1.0); noise = W)
scalar_sol = solve(scalar_prob, SRIW1())
Plots.plot(scalar_sol)

# System of SDEs with Non-Diagonal Noise
function g(du, u, p, t)
    du[1, 1] = 0.3u[1]
    du[1, 2] = 0.6u[1]
    du[1, 3] = 0.9u[1]
    du[1, 4] = 0.12u[1]
    du[2, 1] = 1.2u[2]
    du[2, 2] = 0.2u[2]
    du[2, 3] = 0.3u[2]
    du[2, 4] = 1.8u[2]

    return nothing
end
non_diag_prob = SDEProblem(
    f, g, ones(2), (0.0, 1.0); noise_rate_prototype = zeros(2, 4)
)

non_diag_sol = solve(non_diag_prob)
Plots.plot(non_diag_sol)
