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

# end
