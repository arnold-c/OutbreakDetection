# module NoiseFunctions
#
# export create_noise_arr, create_noise_arr!

# include("ensemble-functions.jl")
# using .EnsembleFunctions

function create_poisson_noise_arr(incarr, noise_spec::NoiseSpecification; seed = 1234)
    noise_arr = zeros(Int64, size(incarr, 1), 1, size(incarr, 3))

    create_poisson_noise_arr!(
        noise_arr, incarr, noise_spec.noise_mean_scaling; seed = seed
    )

    return noise_arr
end

function create_poisson_noise_arr!(
    noise_arr, incarr, noise_mean_scaling; seed = 1234
)
    Random.seed!(seed)
    @inbounds for sim in axes(incarr, 3)
        @views noise_arr[:, 1, sim] .= rand(
            Poisson(noise_mean_scaling * mean(incarr[:, 1, sim])),
            size(incarr, 1),
        )
    end
end

# end
