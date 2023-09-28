# module NoiseFunctions
#
# export create_noise_arr, create_noise_arr!

using DifferentialEquations
using FLoops

# include("ensemble-functions.jl")
# using .EnsembleFunctions

function create_poisson_noise_arr(
    noise_mean, timeparams::SimTimeParameters, nsims
)
    noise_arr = zeros(Int64, timeparams.tlength, 1, nsims)

    create_poisson_noise_arr!(noise_arr, noise_mean, timeparams)

    return noise_arr
end

function create_poisson_noise_arr!(noise_arr, noise_mean, timeparams)
    Random.seed!(1234)
    for sim in axes(noise_arr, 3)
        noise_arr[:, 1, sim] .= rand(Poisson(noise_mean), timeparams.tlength)
    end
end

function create_static_noise_arr(
    init_noise,
    timeparams::SimTimeParameters,
    # dynamicsparams,
    nsims;
    callback = DiscreteCallback(
        sde_condition, sde_affect!; save_positions = (false, false)
    ),
    ode_function = background_ode!,
    noise_function = background_noise!,
)
    # Set noise arr to 3D array (even though not necessary), so it has the same
    # dimensions as the other arrays
    noise_arr = zeros(
        Float64, timeparams.tlength, 1, nsims
    )

    create_static_noise_arr!(
        noise_arr,
        init_noise,
        timeparams;
        # dynamicsparams;
        callback = callback,
        ode_function = ode_function,
        noise_function = noise_function,
    )

    return noise_arr
end

function create_static_noise_arr!(
    noise_arr,
    init_noise::Vector{Float64},
    timeparams::SimTimeParameters;
    # dynamicsparams;
    callback = DiscreteCallback(
        sde_condition, sde_affect!; save_positions = (false, false)
    ),
    ode_function = background_ode!,
    noise_function = background_noise!,
)
    @floop for sim in axes(noise_arr, 3)
        noise_prob = SDEProblem(
            ode_function,
            noise_function,
            init_noise,
            timeparams.tspan,
            # dynamicsparams,
        )

        noise_sol = solve(
            noise_prob, SRIW1(); callback = callback,
            dt = timeparams.tstep,
            adaptive = false,
        )

        noise_arr[:, 1, sim] = @view(noise_sol[1, :])
        # Set first noise incidence to 0 as no new noise in the first time step
        noise_arr[1, 1, sim] = 0.0
    end
    return nothing
end

background_ode!(du, u, p, t) = (du .= 0.0)
background_noise!(du, u, p, t) = (du .= 0.1)

sde_condition(u, t, integrator) = true
function sde_affect!(integrator)
    if integrator.u[1] < 0.0
        integrator.u[1] = -integrator.u[1]
    end
end

# end
