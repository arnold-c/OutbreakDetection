module Helpers
export outdir, getdirpath

using DrWatson: DrWatson

# include("structs.jl")
# using .ODStructs: NoiseSpecification

include("types.jl")
# using .Types: AbstractNoise

outdir(args...) = DrWatson.projectdir("out", args...)

function getdirpath(spec::AbstractNoiseSpecification)
    return reduce(
        joinpath,
        map(
            p -> "$(p)_$(getproperty(spec, p))",
            propertynames(spec),
        )
    )
end

end
